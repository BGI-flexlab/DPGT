package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.ReadLikelihoods;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.GaeaVCFHeaderLines;
import org.bgi.flexlab.gaea.util.GaeaVCFConstants;
import org.bgi.flexlab.gaea.util.Utils;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;

public final class StrandBiasBySample extends GenotypeAnnotation {

    @Override
    public void annotate(final ChromosomeInformationShare ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        if ( likelihoods == null || !g.isCalled() ) {
            return;
        }

        final int[][] table = FisherStrand.getContingencyTable(likelihoods, vc, 0, Arrays.asList(g.getSampleName()));

        gb.attribute(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY, getContingencyArray(table));
    }

    //For now this is only for 2x2 contingency tables
    private static final int ARRAY_DIM = 2;

    /**
     * Helper function to turn the FisherStrand 2x2 table into the SB annotation array
     * @param table the 2x2 table used by the FisherStrand annotation
     * @return the array used by the per-sample Strand Bias annotation
     */
    @VisibleForTesting
    static List<Integer> getContingencyArray(final int[][] table) {
        if(table.length != ARRAY_DIM || table[0].length != ARRAY_DIM) {
            throw new IllegalArgumentException("Expecting a " + ARRAY_DIM + "x" + ARRAY_DIM + " strand bias table.");
        }

        final List<Integer> list = new ArrayList<>(ARRAY_DIM * ARRAY_DIM);
        list.add(table[0][0]);
        list.add(table[0][1]);
        list.add(table[1][0]);
        list.add(table[1][1]);
        return list;
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GaeaVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GaeaVCFHeaderLines.getFormatLine(getKeyNames().get(0)));
    }
}


