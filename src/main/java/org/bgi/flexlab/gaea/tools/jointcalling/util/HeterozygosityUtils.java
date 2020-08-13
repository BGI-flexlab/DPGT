package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.HashMap;
import java.util.Map;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

public class HeterozygosityUtils {

    final public static int REF_INDEX = 0;
    final public static int HET_INDEX = 1;
    final public static int VAR_INDEX = 2;

    protected int sampleCount = -1;
    private Map<Allele, Double> hetCounts;
    private Map<Allele, Double> alleleCounts;
    boolean returnRounded = false;

    /**
     * Create a new HeterozygosityUtils -- a new class should be instantiated for each VariantContext to store data for that VC
     * @param returnRounded round the likelihoods to return integer numbers of counts (as doubles)
     */
    public HeterozygosityUtils(final boolean returnRounded) {
        this.returnRounded = returnRounded;
    }

    /**
     * Get the genotype counts for A/A, A/B, and B/B where A is the reference and B is any alternate allele
     * @param vc
     * @param genotypes may be subset to just founders if a pedigree file is provided
     * @return may be null, otherwise length-3 double[] representing homRef, het, and homVar counts
     */
    public double[] getGenotypeCountsForRefVsAllAlts(final VariantContext vc, final GenotypesContext genotypes) {
        if (genotypes == null || !vc.isVariant())
            return null;

        final boolean doMultiallelicMapping = !vc.isBiallelic();

        int idxAA = 0, idxAB = 1, idxBB = 2;

        double refCount = 0;
        double hetCount = 0;
        double homCount = 0;

        sampleCount = 0;
        for (final Genotype g : genotypes) {
            if (g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2)  // only work for diploid samples
                sampleCount++;
            else
                continue;

            //Need to round the likelihoods to deal with small numerical deviations due to normalizing
            final double[] normalizedLikelihoodsUnrounded = GvcfMathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
            double[] normalizedLikelihoods = new double[normalizedLikelihoodsUnrounded.length];
            if (returnRounded) {
                for (int i = 0; i < normalizedLikelihoodsUnrounded.length; i++) {
                    normalizedLikelihoods[i] = Math.round(normalizedLikelihoodsUnrounded[i]);
                }
            } else {
                normalizedLikelihoods = normalizedLikelihoodsUnrounded;
            }

            if (doMultiallelicMapping) {
                if (g.isHetNonRef()) {
                    //all likelihoods go to homCount
                    homCount++;
                    continue;
                }

                if (!g.isHomRef()) {
                    //get alternate allele for each sample
                    final Allele a1 = g.getAllele(0);
                    final Allele a2 = g.getAllele(1);
                    final int[] idxVector = vc.getGLIndecesOfAlternateAllele(a2.isNonReference() ? a2 : a1);
                    idxAA = idxVector[0];
                    idxAB = idxVector[1];
                    idxBB = idxVector[2];
                }
            }

            refCount += normalizedLikelihoods[idxAA];
            hetCount += normalizedLikelihoods[idxAB];
            homCount += normalizedLikelihoods[idxBB];
        }
        return new double[]{refCount, hetCount, homCount};
    }

    /**
     * Get the count of heterozygotes in vc for a specific altAllele (both reference and non-reference hets, e.g. 1/2)
     * @param vc
     */
    public void doGenotypeCalculations(final VariantContext vc) {
        final GenotypesContext genotypes = vc.getGenotypes();
        if (genotypes == null || !vc.isVariant())
            return;

        final int numAlleles = vc.getNAlleles();

        sampleCount = 0;
        if (hetCounts == null && alleleCounts == null) {
            hetCounts = new HashMap<>();
            alleleCounts = new HashMap<>();
            for (final Allele a : vc.getAlleles()) {
                if (a.isNonReference())
                    hetCounts.put(a, 0.0);
                alleleCounts.put(a, 0.0);
            }

            int idxAB;

            //for each sample
            for (final Genotype g : genotypes) {
                if (g.isCalled() && g.hasLikelihoods() && g.getPloidy() == 2)  // only work for diploid samples
                    sampleCount++;
                else
                    continue;

                int altIndex = 0;
                for(final Allele a : vc.getAlternateAlleles()) {
                    //for each alt allele index from 1 to N
                    altIndex++;

                        final double[] normalizedLikelihoodsUnrounded = GvcfMathUtils.normalizeFromLog10(g.getLikelihoods().getAsVector());
                        double[] normalizedLikelihoods = new double[normalizedLikelihoodsUnrounded.length];
                        if (returnRounded) {
                            for (int i = 0; i < normalizedLikelihoodsUnrounded.length; i++) {
                                normalizedLikelihoods[i] = Math.round(normalizedLikelihoodsUnrounded[i]);
                            }
                        } else {
                            normalizedLikelihoods = normalizedLikelihoodsUnrounded;
                        }
                        //iterate over the other alleles
                        for (int i = 0; i < numAlleles; i++) {
                            //only add homozygotes to alleleCounts, not hetCounts
                            if (i == altIndex) {
                                final double currentAlleleCounts = alleleCounts.get(a);
                                alleleCounts.put(a, currentAlleleCounts + 2*normalizedLikelihoods[GenotypeLikelihoods.calculatePLindex(altIndex,altIndex)]);
                                continue;
                            }
                            //pull out the heterozygote PL index, ensuring that the first allele index < second allele index
                            idxAB = GenotypeLikelihoods.calculatePLindex(Math.min(i,altIndex),Math.max(i,altIndex));
                            final double aHetCounts = hetCounts.get(a);
                            hetCounts.put(a, aHetCounts + normalizedLikelihoods[idxAB]);
                            final double currentAlleleCounts = alleleCounts.get(a);
                            //these are guaranteed to be hets
                            alleleCounts.put(a, currentAlleleCounts + normalizedLikelihoods[idxAB]);
                            final double refAlleleCounts = alleleCounts.get(vc.getReference());
                            alleleCounts.put(vc.getReference(), refAlleleCounts + normalizedLikelihoods[idxAB]);
                        }
                    //add in ref/ref likelihood
                    final double refAlleleCounts = alleleCounts.get(vc.getReference());
                    alleleCounts.put(vc.getReference(), refAlleleCounts + 2*normalizedLikelihoods[0]);
                }

            }
        }
    }

    /**
     * Get the count of heterozygotes in vc for a specific altAllele (both reference and non-reference hets, e.g. 1/2)
     * @param vc
     * @param altAllele the alternate allele of interest
     * @return number of hets
     */
    public double getHetCount(final VariantContext vc, final Allele altAllele) {
        if (hetCounts == null)
            doGenotypeCalculations(vc);
        return hetCounts.containsKey(altAllele)? hetCounts.get(altAllele) : 0;
    }

    public double getAlleleCount(final VariantContext vc, final Allele allele) {
        if (alleleCounts == null)
            doGenotypeCalculations(vc);
        return alleleCounts.containsKey(allele)? alleleCounts.get(allele) : 0;
    }

    public int getSampleCount() {
        return sampleCount;
    }
}
