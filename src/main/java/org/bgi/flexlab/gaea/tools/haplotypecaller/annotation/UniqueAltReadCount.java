package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

public class UniqueAltReadCount extends GenotypeAnnotation {
    public static final String UNIQUE_ALT_READ_SET_COUNT_KEY = "UNIQ_ALT_READ_COUNT";

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(UNIQUE_ALT_READ_SET_COUNT_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(UNIQUE_ALT_READ_SET_COUNT_KEY, 1, VCFHeaderLineType.Integer,
                "Number of ALT reads with unique start and mate end positions at a variant site"));
    }

    @Override
    public void annotate(final ChromosomeInformationShare ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        if (g.isHomRef()) {
            // skip the normal sample
            return;
        }

        final Allele altAllele = vc.getAlternateAllele(0); // assume single-allelic
        final String tumorSampleName = g.getSampleName();
        Collection<ReadLikelihoods<Allele>.BestAllele> tumorBestAlleles = likelihoods.bestAlleles(tumorSampleName);

        // Build a map from the (Start Position, Fragment Size) tuple to the count of reads with that
        // start position and fragment size
        Map<Object, Long> duplicateReadMap = tumorBestAlleles.stream()
                .filter(ba -> ba.allele.equals(altAllele) && ba.isInformative())
                .map(ba -> new ImmutablePair<>(ba.read.getStart(), ba.read.getInferredInsertSize()))
                .collect(Collectors.groupingBy(x -> x, Collectors.counting()));

        gb.attribute(UNIQUE_ALT_READ_SET_COUNT_KEY, duplicateReadMap.size());
    }
}

