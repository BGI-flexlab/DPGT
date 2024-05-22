package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Uses the Exact calculation of Heng Li
 */
abstract class ExactAFCalculator extends AFCalculator {

    protected static final int HOM_REF_INDEX = 0;  // AA likelihoods are always first


    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @param includeDummy   //TODO: Does anyone have any clue what this is?????
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    public static List<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy) {
        return getGLs(GLs, includeDummy, false);
    }

    /**
     * Unpack GenotypesContext into arraylist of double values
     * @param GLs            Input genotype context
     * @param includeDummy   //TODO: Does anyone have any clue what this is?????
     * @param keepUninformative
     * @return               ArrayList of doubles corresponding to GL vectors
     */
    public static List<double[]> getGLs(final GenotypesContext GLs, final boolean includeDummy, boolean keepUninformative) {
        final List<double[]> genotypeLikelihoods = new ArrayList<>(GLs.size() + 1);
        if (includeDummy) {
            genotypeLikelihoods.add((new double[]{0.0, 0.0, 0.0}));
        }

        // Utils.stream(GLs.iterateInSampleNameOrder())
        //         .filter(Genotype::hasLikelihoods)
        //         .map(gt -> gt.getLikelihoods().getAsVector())
        //         .filter(gls -> keepUninformative || GATKVariantContextUtils.isInformative(gls))
        //         .forEach(genotypeLikelihoods::add);
        
        // expand Utils.stream expression
        List<String> sampleNamesInOrder = GLs.getSampleNamesOrderedByName();
        for (String sampleName: sampleNamesInOrder) {
            Genotype gt = GLs.get(sampleName);
            if (gt.hasLikelihoods()) {
                double [] gls = gt.getLikelihoods().getAsVector();
                if (keepUninformative || GATKVariantContextUtils.isInformative(gls)) {
                    genotypeLikelihoods.add(gls);
                }
            }
        }

        return genotypeLikelihoods;
    }

    public static List<int[]> getPLs(final GenotypesContext genotypesContext, final boolean includeDummy)
    {
        return getPLs(genotypesContext, includeDummy, false);
    }

    public static List<int[]> getPLs(final GenotypesContext genotypesContext,
        final boolean includeDummy, boolean keepUninformative)
    {
        List<int[]> PLs = new ArrayList<>(genotypesContext.size() + 1);
        if (includeDummy) {
            PLs.add((new int[]{0, 0, 0}));
        }

        List<String> sampleNamesInOrder = genotypesContext.getSampleNamesOrderedByName();
        for (String sampleName: sampleNamesInOrder) {
            Genotype gt = genotypesContext.get(sampleName);
            if (gt.hasLikelihoods()) {
                int [] pl = gt.getPL();
                if (keepUninformative || GATKVariantContextUtils.isInformative(pl)) {
                    PLs.add(pl);
                }
            }
        }

        return PLs;
    }
}