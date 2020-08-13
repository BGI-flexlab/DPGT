package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.MendelianViolation;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.Trio;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public final class PossibleDeNovo extends InfoFieldAnnotation {

    public PossibleDeNovo(final Set<Trio> trios, final double minGenotypeQualityP) {
        this.trios = Collections.unmodifiableSet(new LinkedHashSet<>(trios));
        if ( trios.isEmpty() ) {
            PossibleDeNovo.logger.warn("Annotation will not be calculated, must provide a valid PED file (-ped) from the command line.");
        }
        mendelianViolation = new MendelianViolation(minGenotypeQualityP);
    }

    /**
     * This is dummy constructor that will do nothing until https://github.com/broadinstitute/gatk/issues/1468 is fixed
     */
    public PossibleDeNovo(){
        this(Collections.emptySet(), 0);
    }

    private static final Logger logger = LogManager.getLogger(PossibleDeNovo.class);

    private static final int hi_GQ_threshold = 20; //WARNING - If you change this value, update the description in GATKVCFHeaderLines
    private static final int lo_GQ_threshold = 10; //WARNING - If you change this value, update the description in GATKVCFHeaderLines
    private static final double percentOfSamplesCutoff = 0.001; //for many, many samples use 0.1% of samples as allele frequency threshold for de novos
    private static final int flatNumberOfSamplesCutoff = 4;

    private final MendelianViolation mendelianViolation;
    private final Set<Trio> trios;

    @Override
    public List<String> getKeyNames() { return Arrays.asList(
            GaeaVCFConstants.HI_CONF_DENOVO_KEY,
            GaeaVCFConstants.LO_CONF_DENOVO_KEY); }

    private static boolean contextHasTrioLikelihoods(final VariantContext vc, final Trio trio) {
        final String mom = trio.getMaternalID();
        final String dad = trio.getPaternalID();
        final String kid = trio.getChildID();

        return   (!mom.isEmpty() && vc.hasGenotype(mom) && vc.getGenotype(mom).hasLikelihoods())
              && (!dad.isEmpty() && vc.hasGenotype(dad) && vc.getGenotype(dad).hasLikelihoods())
              && (!kid.isEmpty() && vc.hasGenotype(kid) && vc.getGenotype(kid).hasLikelihoods());
    }

	@Override
	public Map<String, Object> annotate(ChromosomeInformationShare ref, VariantContext vc,
			ReadLikelihoods<Allele> likelihoods) {
		Utils.nonNull(vc);
        if (trios.isEmpty()){
            return Collections.emptyMap();
        }
        final List<String> highConfDeNovoChildren = new ArrayList<>();
        final List<String> lowConfDeNovoChildren = new ArrayList<>();
        for (final Trio trio : trios) {
            if (vc.isBiallelic() &&
                PossibleDeNovo.contextHasTrioLikelihoods(vc, trio) &&
                mendelianViolation.isViolation(trio.getMother(), trio.getFather(), trio.getChild(), vc) &&
                mendelianViolation.getParentsRefRefChildHet() > 0) {

                final int childGQ = vc.getGenotype(trio.getChildID()).getGQ();
                final int momGQ   = vc.getGenotype(trio.getMaternalID()).getGQ();
                final int dadGQ   = vc.getGenotype(trio.getPaternalID()).getGQ();

                if (childGQ >= hi_GQ_threshold && momGQ >= hi_GQ_threshold && dadGQ >= hi_GQ_threshold) {
                    highConfDeNovoChildren.add(trio.getChildID());
                } else if (childGQ >= lo_GQ_threshold && momGQ > 0 && dadGQ > 0) {
                    lowConfDeNovoChildren.add(trio.getChildID());
                }
            }
        }

        final double percentNumberOfSamplesCutoff = vc.getNSamples()*percentOfSamplesCutoff;
        final double AFcutoff = Math.max(flatNumberOfSamplesCutoff, percentNumberOfSamplesCutoff);
        final int deNovoAlleleCount = vc.getCalledChrCount(vc.getAlternateAllele(0)); //we assume we're biallelic above so use the first alt

        final Map<String,Object> attributeMap = new LinkedHashMap<>(2);
        if ( !highConfDeNovoChildren.isEmpty()  && deNovoAlleleCount < AFcutoff ) {
            attributeMap.put(GaeaVCFConstants.HI_CONF_DENOVO_KEY, highConfDeNovoChildren);
        }
        if ( !lowConfDeNovoChildren.isEmpty()  && deNovoAlleleCount < AFcutoff ) {
            attributeMap.put(GaeaVCFConstants.LO_CONF_DENOVO_KEY, lowConfDeNovoChildren);
        }
        return attributeMap;
	}

}

