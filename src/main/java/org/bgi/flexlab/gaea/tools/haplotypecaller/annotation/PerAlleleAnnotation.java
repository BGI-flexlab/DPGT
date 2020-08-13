package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.stream.Collectors;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.util.QualityUtils;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;

public abstract class PerAlleleAnnotation extends GenotypeAnnotation {

    /**
     * Calculate annotations for eah allele based on given VariantContext and likelihoods for a given genotype's sample
     * and add the annotations to the GenotypeBuilder.  See parent class docs in {@link GenotypeAnnotation}.
     */
    public void annotate(final ChromosomeInformationShare ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(gb);
        Utils.nonNull(vc);
        if ( g == null || likelihoods == null ) {
            return;
        }

        final Map<Allele, List<Integer>> values = likelihoods.alleles().stream()
                .collect(Collectors.toMap(a -> a, a -> new ArrayList<>()));

        Utils.stream(likelihoods.bestAlleles(g.getSampleName()))
                .filter(ba -> ba.isInformative() && isUsableRead(ba.read))
                .forEach(ba -> getValueForRead(ba.read, vc).ifPresent(v -> values.get(ba.allele).add(v)));

        final int[] statistics = vc.getAlleles().stream().filter(this::includeAllele).mapToInt(a -> aggregate(values.get(a))).toArray();
        gb.attribute(getVcfKey(), statistics);
    }

    private boolean includeAllele(final Allele allele) {
        return allele.isNonReference() || includeRefAllele();
    }

    // this is false by default but implementations may wish to override
    protected boolean includeRefAllele() { return false; }

    private static boolean isUsableRead(final GaeaSamRecord read) {
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(getVcfKey(), includeRefAllele() ? VCFHeaderLineCount.R : VCFHeaderLineCount.A, VCFHeaderLineType.Float, getDescription()));
    }

    @Override
    public List<String> getKeyNames() { return Arrays.asList(getVcfKey()); }

    protected abstract OptionalInt getValueForRead(final GaeaSamRecord read, final VariantContext vc);
    protected abstract int aggregate(final List<Integer> values);
    protected abstract String getVcfKey();
    protected abstract String getDescription();
}

