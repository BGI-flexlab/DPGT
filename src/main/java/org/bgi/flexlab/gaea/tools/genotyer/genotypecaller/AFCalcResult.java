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

package org.bgi.flexlab.gaea.tools.genotyer.genotypecaller;

import htsjdk.variant.variantcontext.Allele;
import org.bgi.flexlab.gaea.util.MathUtils;
import org.bgi.flexlab.gaea.util.Utils;

import java.util.*;

/**
 * Describes the results of the AFCalc
 *
 * Only the bare essentials are represented here, as all AFCalc models must return meaningful results for
 * all of these fields.
 *
 * Note that all of the values -- i.e. priors -- are checked now that they are meaningful, which means
 * that users of this code can rely on the values coming out of these functions.
 */
public class AFCalcResult {
    private final static int AF0 = 0;
    private final static int AF1p = 1;
    private final static int LOG_10_ARRAY_SIZES = 2;

    private final double[] log10LikelihoodsOfAC;
    private final double[] log10PriorsOfAC;
    private final double[] log10PosteriorsOfAC;

    private final Map<Allele, Double> log10pRefByAllele;

    /**
     * The AC values for all ALT alleles at the MLE
     */
    private final int[] alleleCountsOfMLE;

    int nEvaluations = 0;

    /**
     * The list of alleles actually used in computing the AF
     */
    private List<Allele> allelesUsedInGenotyping = null;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     */
    public AFCalcResult(final int[] alleleCountsOfMLE,
                        final int nEvaluations,
                        final List<Allele> allelesUsedInGenotyping,
                        final double[] log10LikelihoodsOfAC,
                        final double[] log10PriorsOfAC,
                        final Map<Allele, Double> log10pRefByAllele) {
       
    	if ( allelesUsedInGenotyping == null || allelesUsedInGenotyping.size() < 1 ) throw new IllegalArgumentException("allelesUsedInGenotyping must be non-null list of at least 1 value " + allelesUsedInGenotyping);
        if ( alleleCountsOfMLE == null ) throw new IllegalArgumentException("alleleCountsOfMLE cannot be null");
        if ( alleleCountsOfMLE.length != allelesUsedInGenotyping.size() - 1) throw new IllegalArgumentException("alleleCountsOfMLE.length " + alleleCountsOfMLE.length + " != allelesUsedInGenotyping.size() " + allelesUsedInGenotyping.size());
        if ( nEvaluations < 0 ) throw new IllegalArgumentException("nEvaluations must be >= 0 but saw " + nEvaluations);
        if ( log10LikelihoodsOfAC.length != 2 ) throw new IllegalArgumentException("log10LikelihoodsOfAC must have length equal 2");
        if ( log10PriorsOfAC.length != 2 ) throw new IllegalArgumentException("log10PriorsOfAC must have length equal 2");
        if ( log10pRefByAllele == null ) throw new IllegalArgumentException("log10pRefByAllele cannot be null");
        if ( log10pRefByAllele.size() != allelesUsedInGenotyping.size() - 1 ) throw new IllegalArgumentException("log10pRefByAllele has the wrong number of elements: log10pRefByAllele " + log10pRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        if ( ! allelesUsedInGenotyping.containsAll(log10pRefByAllele.keySet()) ) throw new IllegalArgumentException("log10pRefByAllele doesn't contain all of the alleles used in genotyping: log10pRefByAllele " + log10pRefByAllele + " but allelesUsedInGenotyping " + allelesUsedInGenotyping);
        if ( ! MathUtils.goodLog10ProbVector(log10LikelihoodsOfAC, LOG_10_ARRAY_SIZES, false) ) throw new IllegalArgumentException("log10LikelihoodsOfAC are bad " + Utils.join(",", log10LikelihoodsOfAC));
        if ( ! MathUtils.goodLog10ProbVector(log10PriorsOfAC, LOG_10_ARRAY_SIZES, true) ) throw new IllegalArgumentException("log10priors are bad " + Utils.join(",", log10PriorsOfAC));

        this.alleleCountsOfMLE = alleleCountsOfMLE;
        this.nEvaluations = nEvaluations;
        this.allelesUsedInGenotyping = allelesUsedInGenotyping;

        this.log10LikelihoodsOfAC = Arrays.copyOf(log10LikelihoodsOfAC, LOG_10_ARRAY_SIZES);
        this.log10PriorsOfAC = Arrays.copyOf(log10PriorsOfAC, LOG_10_ARRAY_SIZES);
        this.log10PosteriorsOfAC = computePosteriors(log10LikelihoodsOfAC, log10PriorsOfAC);
        this.log10pRefByAllele = new HashMap<Allele, Double>(log10pRefByAllele);
    }

    /**
     * Return a new AFCalcResult with a new prior probability
     *
     * @param log10PriorsOfAC
     * @return
     */
    public AFCalcResult withNewPriors(final double[] log10PriorsOfAC) {
        return new AFCalcResult(alleleCountsOfMLE, nEvaluations, allelesUsedInGenotyping, log10LikelihoodsOfAC, log10PriorsOfAC, log10pRefByAllele);
    }

    /**
     * Returns a vector with maxAltAlleles values containing AC values at the MLE
     *
     * The values of the ACs for this call are stored in the getAllelesUsedInGenotyping order,
     * starting from index 0 (i.e., the first alt allele is at 0).  The vector is always
     * maxAltAlleles in length, and so only the first getAllelesUsedInGenotyping.size() - 1 values
     * are meaningful.
     *
     * @return a vector with allele counts, not all of which may be meaningful
     */
    //@Ensures("result != null")
    public int[] getAlleleCountsOfMLE() {
        return alleleCountsOfMLE;
    }

    /**
     * Returns the AC of allele a la #getAlleleCountsOfMLE
     *
     * @param allele the allele whose AC we want to know.  Error if its not in allelesUsedInGenotyping
     * @throws IllegalStateException if allele isn't in allelesUsedInGenotyping
     * @return the AC of allele
     */
    public int getAlleleCountAtMLE(final Allele allele) {
        return getAlleleCountsOfMLE()[altAlleleIndex(allele)];
    }

    /**
     * Returns the number of cycles used to evaluate the pNonRef for this AF calculation
     *
     * @return the number of evaluations required to produce the answer for this AF calculation
     */
    public int getnEvaluations() {
        return nEvaluations;
    }

    /**
     * Get the list of alleles actually used in genotyping.
     *
     * Due to computational / implementation constraints this may be smaller than
     * the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping, the first of which is the reference allele
     */
    //@Ensures({"result != null", "! result.isEmpty()"})
    public List<Allele> getAllelesUsedInGenotyping() {
        return allelesUsedInGenotyping;
    }

    /**
     * Get the log10 normalized -- across all ACs -- posterior probability of AC == 0 for all alleles
     *
     * @return
     */
    //@Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PosteriorOfAFEq0() {
        return log10PosteriorsOfAC[AF0];
    }

    /**
     * Get the log10 normalized -- across all ACs -- posterior probability of AC > 0 for any alleles
     *
     * @return
     */
    //@Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PosteriorOfAFGT0() {
        return log10PosteriorsOfAC[AF1p];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- likelihood of AC == 0 for all alleles
     *
     * @return
     */
    //@Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10LikelihoodOfAFEq0() {
        return log10LikelihoodsOfAC[AF0];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- likelihood of AC > 0 for any alleles
     *
     * @return
     */
    //@Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10LikelihoodOfAFGT0() {
        return log10LikelihoodsOfAC[AF1p];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- prior probability of AC == 0 for all alleles
     *
     * @return
     */
    //@Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PriorOfAFEq0() {
        return log10PriorsOfAC[AF0];
    }

    /**
     * Get the log10 unnormalized -- across all ACs -- prior probability of AC > 0
     *
     * @return
     */
   // @Ensures({"MathUtils.goodLog10Probability(result)"})
    public double getLog10PriorOfAFGT0() {
        return log10PriorsOfAC[AF1p];
    }

    @Override
    public String toString() {
        final List<String> byAllele = new LinkedList<String>();
        for ( final Allele a : getAllelesUsedInGenotyping() )
            if ( a.isNonReference() ) byAllele.add(String.format("%s => MLE %d / posterior %.2f", a, getAlleleCountAtMLE(a), getLog10PosteriorOfAFEq0ForAllele(a)));
        return String.format("AFCalc%n\t\tlog10PosteriorOfAFGT0=%.2f%n\t\t%s", getLog10LikelihoodOfAFGT0(), Utils.join("\n\t\t", byAllele));
    }

    /**
     * Are we sufficiently confidence in being non-ref that the site is considered polymorphic?
     *
     * We are non-ref if the probability of being non-ref > the emit confidence (often an argument).
     * Suppose posterior AF > 0 is log10: -5 => 10^-5
     * And that log10minPNonRef is -3.
     * We are considered polymorphic since 10^-5 < 10^-3 => -5 < -3
     *
     * Note that log10minPNonRef is really the minimum confidence, scaled as an error rate, so
     * if you want to be 99% confidence, then log10PNonRef should be log10(0.01) = -2.
     *
     * @param log10minPNonRef the log10 scaled min pr of being non-ref to be considered polymorphic
     *
     * @return true if there's enough confidence (relative to log10minPNonRef) to reject AF == 0
     */
    //@Requires("MathUtils.goodLog10Probability(log10minPNonRef)")
    public boolean isPolymorphic(final Allele allele, final double log10minPNonRef) {
        return getLog10PosteriorOfAFEq0ForAllele(allele) < log10minPNonRef;
    }

    /**
     * Same as #isPolymorphic but takes a phred-scaled quality score as input
     */
    public boolean isPolymorphicPhredScaledQual(final Allele allele, final double minPNonRefPhredScaledQual) {
        if ( minPNonRefPhredScaledQual < 0 ) throw new IllegalArgumentException("phredScaledQual " + minPNonRefPhredScaledQual + " < 0 ");
        final double log10Threshold = minPNonRefPhredScaledQual / -10;
        return isPolymorphic(allele, log10Threshold);
    }

    /**
     * Are any of the alleles polymorphic w.r.t. #isPolymorphic?
     *
     * @param log10minPNonRef the confidence threshold, in log10 space
     * @return true if any are poly, false otherwise
     */
    public boolean anyPolymorphic(final double log10minPNonRef) {
        for ( final Allele a : getAllelesUsedInGenotyping() )
            if ( a.isNonReference() && isPolymorphic(a, log10minPNonRef) )
                return true;
        return false;
    }

    /**
     * Returns the log10 probability that allele is not segregating
     *
     * Note that this function is p not segregating so that we can store
     * internally the log10 value of AF == 0, which grows very quickly
     * negative and yet has sufficient resolution for high confidence tests.
     * For example, if log10pRef == -100, not an unreasonably high number,
     * if we tried to store log10pNonRef we'd be looking at 1 - 10^-100, which
     * quickly underflows to 1.  So the logic here is backward from what
     * you really want (the p of segregating) but we do that for numerical
     * reasons
     *
     * Unlike the sites-level annotation, this calculation is specific to allele, and can be
     * used to separately determine how much evidence there is that allele is independently
     * segregating as opposed to the site being polymorphic with any allele.  In the bi-allelic
     * case these are obviously the same but for multiple alt alleles there can be lots of
     * evidence for one allele but not so much for any other allele
     *
     * @param allele the allele we're interested in, must be in getAllelesUsedInGenotyping
     * @return the log10 probability that allele is not segregating at this site
     */
    //@Ensures("MathUtils.goodLog10Probability(result)")
    public double getLog10PosteriorOfAFEq0ForAllele(final Allele allele) {
        final Double log10pNonRef = log10pRefByAllele.get(allele);
        if ( log10pNonRef == null ) throw new IllegalArgumentException("Unknown allele " + allele);
        return log10pNonRef;
    }

    /**
     * Returns the log10 normalized posteriors given the log10 likelihoods and priors
     *
     * @param log10LikelihoodsOfAC
     * @param log10PriorsOfAC
     *
     * @return freshly allocated log10 normalized posteriors vector
     */
    //@Requires("log10LikelihoodsOfAC.length == log10PriorsOfAC.length")
    //@Ensures("MathUtils.goodLog10ProbVector(result, LOG_10_ARRAY_SIZES, true)")
    private static double[] computePosteriors(final double[] log10LikelihoodsOfAC, final double[] log10PriorsOfAC) {
        final double[] log10UnnormalizedPosteriors = new double[log10LikelihoodsOfAC.length];
        for ( int i = 0; i < log10LikelihoodsOfAC.length; i++ )
            log10UnnormalizedPosteriors[i] = log10LikelihoodsOfAC[i] + log10PriorsOfAC[i];
        return MathUtils.normalizeFromLog10(log10UnnormalizedPosteriors, true, false);
    }

    /**
     * Computes the offset into linear vectors indexed by alt allele for allele
     *
     * Things like our MLE allele count vector are indexed by alt allele index, with
     * the first alt allele being 0, the second 1, etc.  This function computes the index
     * associated with allele.
     *
     * @param allele the allele whose alt index we'd like to know
     * @throws IllegalArgumentException if allele isn't in allelesUsedInGenotyping
     * @return an index value greater than 0 suitable for indexing into the MLE and other alt allele indexed arrays
     */
    //@Requires("allele != null")
    //@Ensures({"result >= 0", "result < allelesUsedInGenotyping.size() - 1"})
    private int altAlleleIndex(final Allele allele) {
        if ( allele.isReference() ) throw new IllegalArgumentException("Cannot get the alt allele index for reference allele " + allele);
        final int index = allelesUsedInGenotyping.indexOf(allele);
        if ( index == -1 )
            throw new IllegalArgumentException("could not find allele " + allele + " in " + allelesUsedInGenotyping);
        else
            return index - 1;
    }
}