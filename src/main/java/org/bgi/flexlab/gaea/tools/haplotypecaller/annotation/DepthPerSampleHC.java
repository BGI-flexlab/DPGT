package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public final class DepthPerSampleHC extends GenotypeAnnotation implements StandardHCAnnotation {

    @Override
    public void annotate( final ChromosomeInformationShare ref,
                          final VariantContext vc,
                          final Genotype g,
                          final GenotypeBuilder gb,
                          final ReadLikelihoods<Allele> likelihoods ) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        if ( likelihoods == null || !g.isCalled() ) {
            return;
        }

        // check that there are reads
        final String sample = g.getSampleName();
        if (likelihoods.sampleReadCount(likelihoods.indexOfSample(sample)) == 0) {
            gb.DP(0);
            return;
        }

        final Set<Allele> alleles = new LinkedHashSet<>(vc.getAlleles());

        // make sure that there's a meaningful relationship between the alleles in the likelihoods and our VariantContext
        if ( !likelihoods.alleles().containsAll(alleles) ) {
            return;
        }

        // the depth for the HC is the sum of the informative alleles at this site.  It's not perfect (as we cannot
        // differentiate between reads that align over the event but aren't informative vs. those that aren't even
        // close) but it's a pretty good proxy and it matches with the AD field (i.e., sum(AD) = DP).
        final Map<Allele, List<Allele>> alleleSubset = alleles.stream().collect(Collectors.toMap(a -> a, a -> Arrays.asList(a)));
        final ReadLikelihoods<Allele> subsettedLikelihoods = likelihoods.marginalize(alleleSubset);
        final int depth = (int) subsettedLikelihoods.bestAlleles(sample).stream().filter(ba -> ba.isInformative()).count();
        gb.DP(depth);
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(VCFConstants.DEPTH_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
    }
}

