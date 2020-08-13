package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
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

public final class DepthPerAlleleBySample extends GenotypeAnnotation implements StandardAnnotation, StandardMutectAnnotation {

    @Override
    public void annotate(final ChromosomeInformationShare ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb, "gb is null");
        Utils.nonNull(vc, "vc is null");

        if ( g == null || !g.isCalled() || likelihoods == null) {
            return;
        }
        final Set<Allele> alleles = new LinkedHashSet<>(vc.getAlleles());

        // make sure that there's a meaningful relationship between the alleles in the likelihoods and our VariantContext
        Utils.validateArg(likelihoods.alleles().containsAll(alleles), () -> "VC alleles " + alleles + " not a  subset of ReadLikelihoods alleles " + likelihoods.alleles());

        final Map<Allele, Integer> alleleCounts = new LinkedHashMap<>();
        for ( final Allele allele : vc.getAlleles() ) {
            alleleCounts.put(allele, 0);
        }
        final Map<Allele, List<Allele>> alleleSubset = alleles.stream().collect(Collectors.toMap(a -> a, Arrays::asList));
        final ReadLikelihoods<Allele> subsettedLikelihoods = likelihoods.marginalize(alleleSubset);
        subsettedLikelihoods.bestAlleles(g.getSampleName()).stream()
                .filter(ba -> ba.isInformative())
                .forEach(ba -> alleleCounts.compute(ba.allele, (allele,prevCount) -> prevCount + 1));

        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference()); //first one in AD is always ref
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++) {
            counts[i + 1] = alleleCounts.get(vc.getAlternateAllele(i));
        }

        gb.AD(counts);
    }

    @Override
    public List<String> getKeyNames() { return Collections.singletonList(VCFConstants.GENOTYPE_ALLELE_DEPTHS); }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(VCFStandardHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}
