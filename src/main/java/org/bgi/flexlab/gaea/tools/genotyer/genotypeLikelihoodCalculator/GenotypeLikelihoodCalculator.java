/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator;

import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.mapreduce.genotyper.GenotyperOptions;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by zhangyong on 2016/12/21.
 * Genotype likelihood calculators, most source code came from GATK UnifiedGenotyperEngine Class
 */
public abstract class GenotypeLikelihoodCalculator {
    /**Class
     * Calculator Name
     */
    public static String GTMSTRING = "GenotypeLikelihoodCalculator";

    /**
     * models for genotype calculator
     */
    public enum Model {
        SNP,
        INDEL,
        SOAPSNP,
        SamtoolsINDEL,
        BOTH
    }

    /**
     * models to be used
     */
    public static List<Model> modelsToUse = null;

    /**
     * name -> genotype likelihood calculator objects
     */
    public static Map<String, GenotypeLikelihoodCalculator> glcm = new HashMap<>();

    /**
     * get genotype likelihood calculators models to be used from options
     * @param options options
     * @return genotype likelihood model list
     */
    public static List<Model> getCalculators(GenotyperOptions options) {
        if(options.getGtlcalculators() == Model.BOTH) {
            modelsToUse = new ArrayList<>(2);
            modelsToUse.add(Model.SNP);
            modelsToUse.add(Model.INDEL);
        } else {
            modelsToUse = new ArrayList<>(1);
            modelsToUse.add(options.getGtlcalculators());
        }
            return modelsToUse;
    }

    /**
     * get map -> genotype likelihood calculator objects
     * @param options options
     * @return genotype likelihood calculator Classes
     */
    public static Map<String, GenotypeLikelihoodCalculator> getGenotypeLikelihoodsCalculatorObject(GenotyperOptions options) {
        List<Class<? extends GenotypeLikelihoodCalculator>> glmClasses = new ArrayList<>();
        glmClasses.add(SNPGenotypeLikelihoodCalculator.class);
        glmClasses.add(INDELGenotypeLikelihoodCalculator.class);
        //glmClasses.add(SOAPSNPGenotypeLikelihoodCalculator.class);
        for (int i = 0; i < glmClasses.size(); i++) {
            final Class<? extends GenotypeLikelihoodCalculator> glmClass = glmClasses.get(i);
            final String key = glmClass.getSimpleName().replaceAll(GTMSTRING, "").toUpperCase();
            try {
                final Object args[] = new Object[] { options };
                final Constructor c = glmClass.getDeclaredConstructor(GenotyperOptions.class);
                glcm.put(key, (GenotypeLikelihoodCalculator) c.newInstance(args));
                System.err.println("key:" + key + "\tclass:" + c.getName());
            } catch (Exception e) {
                throw new UserException("The likelihoods model provided for the -glm argument (" + options.getGtlcalculators() +
                        ") is not a valid option: " + e.getClass() + ":\n" + e.getMessage());
            }
        }

        return glcm;
    }

    /**
     * constructor
     * @param options options
     */
    public GenotypeLikelihoodCalculator(GenotyperOptions options) {
        if(options == null) {
            throw new RuntimeException("options can not be null.");
        }
    }

    /**
     * Alleles number for Cache
     */
    private final static int NUM_LIKELIHOODS_CACHE_N_ALLELES = 5;

    /**
     * ploidy number for Cache
     */
    private final static int NUM_LIKELIHOODS_CACHE_PLOIDY = 10;

    /**
     * the number of likelihoods or genotype
     * caching numAlleles up to 5 and ploidy up to 10
     */
    private final static int[][] numLikelihoodCache = new int[NUM_LIKELIHOODS_CACHE_N_ALLELES][NUM_LIKELIHOODS_CACHE_PLOIDY];

    /**
     * max phred-likelihood (PL)
     */
    public final static int MAX_PL = Short.MAX_VALUE;

    /**
     * The maximum number of alleles that we can represent as genotype likelihoods
     */
    public final static int MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED = 50;

    /**
    * a cache of the PL index to the 2 alleles it represents over all possible numbers of alternate alleles
    */
    private final static GenotypeLikelihoodsAllelePair[] PLIndexToAlleleIndex = calculatePLcache(MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);

    /**
     * Compute how many likelihood elements are associated with the given number of alleles
     * Equivalent to asking in how many ways N non-negative integers can add up to P is S(N,P)
     * where P = ploidy (number of chromosomes) and N = total # of alleles.
     * Each chromosome can be in one single state (0,...,N-1) and there are P of them.
     * Naive solution would be to store N*P likelihoods, but this is not necessary because we can't distinguish chromosome states, but rather
     * only total number of alt allele counts in all chromosomes.
     *
     * For example, S(3,2) = 6: For alleles A,B,C, on a diploid organism we have six possible genotypes:
     * AA,AB,BB,AC,BC,CC.
     * Another way of expressing is with vector (#of A alleles, # of B alleles, # of C alleles)
     * which is then, for ordering above, (2,0,0), (1,1,0), (0,2,0), (1,1,0), (0,1,1), (0,0,2)
     * In general, for P=2 (regular biallelic), then S(N,2) = N*(N+1)/2
     *
     * Note this method caches the value for most common num Allele / ploidy combinations for efficiency
     *
     * Recursive implementation:
     *   S(N,P) = sum_{k=0}^P S(N-1,P-k)
     *  because if we have N integers, we can condition 1 integer to be = k, and then N-1 integers have to sum to P-K
     * With initial conditions
     *   S(N,1) = N  (only way to have N integers add up to 1 is all-zeros except one element with a one. There are N of these vectors)
     *   S(1,P) = 1 (only way to have 1 integer add to P is with that integer P itself).
     *
     *   @param  numAlleles      Number of alleles (including ref)
     *   @param  ploidy          Ploidy, or number of chromosomes in set
     *   @return    Number of likelihood elements we need to hold.
     */
    public static int numLikelihoods(final int numAlleles, final int ploidy) {
        if ( numAlleles < NUM_LIKELIHOODS_CACHE_N_ALLELES
                && ploidy < NUM_LIKELIHOODS_CACHE_PLOIDY )
            return numLikelihoodCache[numAlleles][ploidy];
        else {
            // have to calculate on the fly
            return calcNumLikelihoods(numAlleles, ploidy);
        }
    }

    /**
     * Actually does the computation in @see #numLikelihoods
     *
     * @param numAlleles number of alleles
     * @param ploidy number of ploidy
     * @return number of likelihoods
     */
    private static final int calcNumLikelihoods(final int numAlleles, final int ploidy) {
        if (numAlleles == 1)
            return 1;
        else if (ploidy == 1)
            return numAlleles;
        else {
            int acc =0;
            for (int k=0; k <= ploidy; k++ )
                acc += calcNumLikelihoods(numAlleles - 1, ploidy - k);
            return acc;
        }
    }

    /**
     * get the allele index pair for the given PL
     *
     * @param PLindex   the PL index
     * @return the allele index pair
     */
    public static GenotypeLikelihoodsAllelePair getAllelePair(final int PLindex) {
        // make sure that we've cached enough data
        if ( PLindex >= PLIndexToAlleleIndex.length )
            throw new UserException("Gaea limitation: cannot genotype more than " + MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED + " alleles");

        return PLIndexToAlleleIndex[PLindex];
    }

    /**
     * As per the VCF spec: "the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
     * In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
     * Assumes that allele1Index < allele2Index
     */
    public static int calculatePLindex(final int allele1Index, final int allele2Index) {
        return (allele2Index * (allele2Index + 1) / 2) + allele1Index;
    }

    /**
     *
     * @param altAlleles alt alleles number
     * @return allele pair of all possible
     */
    private static GenotypeLikelihoodsAllelePair[] calculatePLcache(final int altAlleles) {
        final int numLikelihoods = numLikelihoods(1 + altAlleles, 2);
        final GenotypeLikelihoodsAllelePair[] cache = new GenotypeLikelihoodsAllelePair[numLikelihoods];

        // for all possible combinations of 2 alleles
        for ( int allele1 = 0; allele1 <= altAlleles; allele1++ ) {
            for ( int allele2 = allele1; allele2 <= altAlleles; allele2++ ) {
                cache[calculatePLindex(allele1, allele2)] = new GenotypeLikelihoodsAllelePair(allele1, allele2);
            }
        }

        // a bit of sanity checking
        for ( int i = 0; i < cache.length; i++ ) {
            if ( cache[i] == null )
                throw new UserException("BUG: cache entry " + i + " is unexpected null");
        }

        return cache;
    }

    /**
     *  Genotype Allele Pair Class
     */
    public static class GenotypeLikelihoodsAllelePair {
        public final int alleleIndex1, alleleIndex2;

        public GenotypeLikelihoodsAllelePair(final int alleleIndex1, final int alleleIndex2) {
            this.alleleIndex1 = alleleIndex1;
            this.alleleIndex2 = alleleIndex2;
        }
    }

    /**
     * abstract method for genotype likelihood calculation
     * @param mpileup multi-sample pileups
     * @param reference reference
     * @param options options
     * @return variantContext with Likelihoods
     */
    public abstract VariantContext genotypeLikelihoodCalculate(Mpileup mpileup, ChromosomeInformationShare reference, GenotyperOptions options, GenomeLocationParser locationParser, Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap);
}
