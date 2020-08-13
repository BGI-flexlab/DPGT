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

import htsjdk.variant.variantcontext.*;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.genotyer.VariantCallingEngine;
import org.bgi.flexlab.gaea.tools.mapreduce.genotyper.GenotyperOptions;
import org.bgi.flexlab.gaea.util.BaseUtils;
import org.bgi.flexlab.gaea.util.DiploidGenotype;
import org.bgi.flexlab.gaea.util.MathUtils;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static org.apache.commons.math3.util.FastMath.log10;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * Created by zhangyong on 2016/12/21.
 */
public class SNPGenotypeLikelihoodCalculator extends GenotypeLikelihoodCalculator {
    /**
     * Constant static data: log3
     */
    protected final static double log10_3 = log10(3.0);

    /**
     * Constant static data: default pcr error rate
     */
    public final static double DEFAULT_PCR_ERROR_RATE = 1e-4;

    /**
     * Constant static data: ploidy, 2 as default
     */
    protected final static int FIXED_PLOIDY = 2;

    /**
     * Constant static data: max ploid
     */
    protected final static int MAX_PLOIDY = FIXED_PLOIDY + 1;

    /**
     * Constant static data: ploidy adjustment for likelihood calculation
     */
    protected final static double ploidyAdjustment = log10(FIXED_PLOIDY);

    /**
     * base quality calculation cache
     * FIXME::only consider diploid and without overlapping paired reads
     */
    public static GenotypeData[][] CACHE = new GenotypeData[BaseUtils.BASES.length][QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE+1];

    /**
     * one time calculation result about PCR error
     */
    protected static double log10_PCR_error_3 = log10(DEFAULT_PCR_ERROR_RATE) - log10_3;
    protected static double log10_1_minus_PCR_error = log10(1.0 - DEFAULT_PCR_ERROR_RATE);


    /**
     * get Cached likelihoods
     * @param base base Index
     * @param baseQuality base quality
     * @return likelihoods of this base to each genotype
     */
    public GenotypeData getCACHELikelihoods(byte base, byte baseQuality) {
        return CACHE[base][baseQuality];
    }

    /**
     * is likelihoods of this base to each genotype in Cache
     * @param base base Index
     * @param baseQuality base quality
     * @return is in cache
     */
    public boolean hasCACHE(byte base, byte baseQuality) {
        return getCACHELikelihoods(base, baseQuality) != null;
    }

    /**
     * set cache
     * @param base base Index
     * @param baseQuality base quality
     * @param genotypeData likelihoods of this base to each genotype
     */
    public void setCACHE(byte base, byte baseQuality, GenotypeData genotypeData) {
        CACHE[base][baseQuality] = genotypeData;
    }

    /**
     * sum likelihoods of each allele
     */
    private final double[] likelihoodSums = new double[4];

    /**
     * constructor
     * @param options
     */
    public SNPGenotypeLikelihoodCalculator(GenotyperOptions options) {
        super(options);
        init(options.getPcr_error());
    }


    /**
     * constructor
     * @param options options
     * @param pcrErrorRate pcr error rate
     */
    public SNPGenotypeLikelihoodCalculator(GenotyperOptions options, double pcrErrorRate) {
        super(options);
        init(pcrErrorRate);
    }

    public void init(double pcrErrorRate) {
        if(pcrErrorRate == 0) {
            pcrErrorRate = DEFAULT_PCR_ERROR_RATE;
        }
        log10_PCR_error_3 = log10(pcrErrorRate) - log10_3;
        log10_1_minus_PCR_error = log10(1.0 - pcrErrorRate);
    }

    /**
     * calculate genotype likelihoods for SNP
     * @param mpileup multi-sample pileups
     * @param reference reference
     * @param options options
     * @return variant context with likelihoods
     */
    public VariantContext genotypeLikelihoodCalculate(Mpileup mpileup, ChromosomeInformationShare reference, GenotyperOptions options, GenomeLocationParser locationParser, final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        final byte refBase = (byte)reference.getBase(mpileup.getPosition());
        final int indexOfRefBase = BaseUtils.simpleBaseToBaseIndex(refBase);
        // handle non-standard reference bases
        if ( indexOfRefBase == -1 )
            return null;
        final Allele refAllele = Allele.create(refBase, true);

        // calculate the GLs
        ArrayList<SampleGenotypeData> gls = new ArrayList<>(mpileup.getSize());
        Map<String, Pileup> pileups = mpileup.getCurrentPosPileup();
        int position = mpileup.getPosition();
        if (pileups != null) {
            //calculate the genotype likelihood
            //System.err.println("ref:" + (char) refBase + "\tat:" + position);
            for(String sample : pileups.keySet()) {
                Pileup pileup = pileups.get(sample);
                //System.err.println(sample);
                //depth too low to calculate genotype likelihood
                if (pileup.getFilteredPileup().size() < options.getMinDepth() || pileup.getDeletionRate() > 0.2) {
                    continue;
                }

                //calculation genotype likelihoods
                SampleGenotypeData sampleGenotypeData = getGenotypeLikelihood(pileup, options.isCapBaseQualsAtMappingQual(), options.getMinBaseQuality());
                if(sampleGenotypeData.getDepth() > options.getMinDepth()) {
                    sampleGenotypeData.setName(sample);
                    gls.add(sampleGenotypeData);
                }
                //System.err.println("genotype likelihood result:" + sampleGenotypeData.toString());
            }
        }

        //build variant context
        List<Allele> alleles = new ArrayList<>();
        alleles.add(refAllele);
        VariantContextBuilder builder = new VariantContextBuilder("GaeaCall", reference.getChromosomeName(),
                position + 1, position + 1, alleles);

        //System.err.println("determine alt alleles:");
        alleles.addAll(determineAlternateAlleles(refBase, gls));
        // if there are no non-ref alleles...
        if ( alleles.size() == 1 ) {
            if ( options.getOutputMode() == VariantCallingEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                return builder.make();

            // otherwise, choose any alternate allele (it doesn't really matter)
            alleles.add(Allele.create(BaseUtils.baseIndexToSimpleBase(indexOfRefBase == 0 ? 1 : 0)));
        }

        // create the alternate alleles and the allele ordering (the ordering is crucial for the GLs)
        final int numAlleles = alleles.size();
        final int numAltAlleles = numAlleles - 1;

        final int[] alleleOrdering = new int[numAlleles];
        int alleleOrderingIndex = 0;
        int numLikelihoods = 0;
        for ( Allele allele : alleles ) {
            alleleOrdering[alleleOrderingIndex++] = BaseUtils.simpleBaseToBaseIndex(allele.getBases()[0]);
            numLikelihoods += alleleOrderingIndex;
        }
        builder.alleles(alleles);

        // create the PL ordering to use based on the allele ordering.
        final int[] PLordering = new int[numLikelihoods];
        for ( int i = 0; i <= numAltAlleles; i++ ) {
            for ( int j = i; j <= numAltAlleles; j++ ) {
                // As per the VCF spec: "the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
                // In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
                PLordering[(j * (j+1) / 2) + i] = DiploidGenotype.createDiploidGenotype(alleleOrdering[i], alleleOrdering[j]).ordinal();
            }
        }

        // create the genotypes; no-call everyone for now
        final GenotypesContext genotypes = GenotypesContext.create();

        for ( SampleGenotypeData sampleData : gls ) {
            final double[] allLikelihoods = sampleData.getLog10Likelihoods();
            final double[] myLikelihoods = new double[numLikelihoods];

            for ( int i = 0; i < numLikelihoods; i++ )
                myLikelihoods[i] = allLikelihoods[PLordering[i]];

            // normalize in log space so that max element is zero.
            final GenotypeBuilder gb = new GenotypeBuilder(sampleData.getName());
            final double[] genotypeLikelihoods = MathUtils.normalizeFromLog10(myLikelihoods, false, true);
            gb.PL(genotypeLikelihoods);
            gb.DP(sampleData.getDepth());
            genotypes.add(gb.make());
        }

        return builder.genotypes(genotypes).make();
    }

    /**
     * determines the alleles to use
     * @param ref ref base
     * @param sampleDataList sample genotype data
     * @return allele list
     */
    protected List<Allele> determineAlternateAlleles(final byte ref, final List<SampleGenotypeData> sampleDataList) {

        final int baseIndexOfRef = BaseUtils.simpleBaseToBaseIndex(ref);
        final int PLindexOfRef = DiploidGenotype.createDiploidGenotype(ref, ref).ordinal();
        for ( int i = 0; i < 4; i++ )
            likelihoodSums[i] = 0.0;

        // based on the GLs, find the alternate alleles with enough probability
        for ( SampleGenotypeData sampleData : sampleDataList ) {
            //System.err.println(sampleData.toString());
            final double[] likelihoods = sampleData.getLog10Likelihoods();
            final int PLindexOfBestGL = MathUtils.maxElementIndex(likelihoods);
            //System.err.println("ref index:" + PLindexOfRef + "\t best genotype index:" + PLindexOfBestGL);
            if ( PLindexOfBestGL != PLindexOfRef ) {
                GenotypeLikelihoodsAllelePair alleles = getAllelePair(PLindexOfBestGL);
               // System.err.println("best genotype allele1:" + (char) BaseUtils.baseIndexToSimpleBase(alleles.alleleIndex1) +
               //         "\tbest genotype allele2:" + (char) BaseUtils.baseIndexToSimpleBase(alleles.alleleIndex2));
                if ( alleles.alleleIndex1 != baseIndexOfRef )
                    likelihoodSums[alleles.alleleIndex1] += likelihoods[PLindexOfBestGL] - likelihoods[PLindexOfRef];
                // don't double-count it
                if ( alleles.alleleIndex2 != baseIndexOfRef && alleles.alleleIndex2 != alleles.alleleIndex1 )
                    likelihoodSums[alleles.alleleIndex2] += likelihoods[PLindexOfBestGL] - likelihoods[PLindexOfRef];
            }
        }

        final List<Allele> allelesToUse = new ArrayList<>();
        for ( int i = 0; i < 4; i++ ) {
            if ( likelihoodSums[i] > 0.0 ) {
                //System.err.println("allele:" + (char) BaseUtils.baseIndexToSimpleBase(i) + "\tlikelihood sum:" + likelihoodSums[i]);
                allelesToUse.add(Allele.create(BaseUtils.baseIndexToSimpleBase(i), false));
            }
        }

        return allelesToUse;
    }

    /**
     * calculate the genotype likelihood for one sample
     * @param pileup pileup
     * @param isCapBaseQualsAtMappingQual options
     * @return sample genotype likelihoods data
     */
    private SampleGenotypeData getGenotypeLikelihood(Pileup pileup, boolean isCapBaseQualsAtMappingQual, int minBaseQuality) {
        int goodBaseCount = 0;
        SampleGenotypeData sampleGenotypeData = new SampleGenotypeData();
       // System.err.println("depth:" + pileup.getTotalPileup().size());
        for(PileupReadInfo readInfo : pileup.getTotalPileup()) {
            if(readInfo.isDeletionBase())
                continue;

            byte base = readInfo.getBinaryBase();
            byte quality = readInfo.getBaseQuality();

            //System.err.println("base before:" + (char) BaseUtils.baseIndexToSimpleBase(base) + "\tquality:" + (char) (quality + 33));

            if(isCapBaseQualsAtMappingQual && quality > readInfo.getMappingQuality()) {
                quality = (byte)readInfo.getMappingQuality();
            }
            if(quality < minBaseQuality)
                quality = 0;

            if(base < 0 || base > 3) {
                //System.err.println("bad base" + readInfo.getBase());
                continue;
            }

            //System.err.println("base after:" + (char) BaseUtils.baseIndexToSimpleBase(base) + "\tquality:" + quality);

            goodBaseCount++;
            if(hasCACHE(base, quality)) {
                sampleGenotypeData.add(getCACHELikelihoods(base, quality));
            } else {
                GenotypeData genotypeDataTmp = errorModel(base, quality);
                sampleGenotypeData.add(genotypeDataTmp);
                setCACHE(base, quality, genotypeDataTmp);
            }
        }

        sampleGenotypeData.setDepth(goodBaseCount);

        return sampleGenotypeData;
    }

    /**
     * calculate genotype likelihood
     * FIXME::only consider diploid and without overlapping paired reads
     * @param observedBinaryBase base
     * @param baseQuality quality
     * @return genotype likelihood
     */
    private GenotypeData errorModel(byte observedBinaryBase, byte baseQuality) {
        byte observedBase = BaseUtils.baseIndexToSimpleBase(observedBinaryBase);
        GenotypeData genotypeData = new GenotypeData();
        double[] log10FourBaseLikelihoods = baseErrorModel(observedBase, baseQuality);
        for ( DiploidGenotype g : DiploidGenotype.values() ) {

            // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
            double p_base = 0.0;
            p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - ploidyAdjustment);
            p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - ploidyAdjustment);
            final double likelihood = log10(p_base);
            //System.err.println("genotype:" + g.name() + "\tlikelihood:" + p_base);
            genotypeData.log10Likelihoods[g.ordinal()] += likelihood;
        }
        return genotypeData;
    }

    /**
     * base error likelihood
     * @param observedBase base
     * @param baseQuality quality
     * @return base error likelihood
     */
    private double[] baseErrorModel(byte observedBase, byte baseQuality) {
        double[] log10FourBaseLikelihoods = GenotypeData.baseZeros.clone();
        for (byte trueBase : BaseUtils.BASES) {
            double likelihood = 0.0;

            //System.err.println("true base:" + (char) trueBase);
            for (byte fragmentBase : BaseUtils.BASES) {
                double log10FragmentLikelihood = (trueBase == fragmentBase ? log10_1_minus_PCR_error : log10_PCR_error_3);
                //System.err.print("\tfragment base:" + (char)fragmentBase + "\tpcr error prob:" + log10FragmentLikelihood);
                if(baseQuality != 0) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase, fragmentBase, baseQuality);
                }
                //System.err.println("\tquality prob:" + log10PofObservingBaseGivenChromosome(observedBase, fragmentBase, baseQuality) + "\tsum:" + log10FragmentLikelihood);
                likelihood += pow(10, log10FragmentLikelihood);
            }
            //System.err.println("\tfinal prob:" + log10(likelihood));
            log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(trueBase)] = log10(likelihood);
        }
        return log10FourBaseLikelihoods;
    }

    /**
     * calculate quality error prob in log10
     * @param observedBase observed base
     * @param fragmentBase fragment base
     * @param baseQuality quality
     * @return quality error prob in log10
     */
    private double log10PofObservingBaseGivenChromosome(byte observedBase, byte fragmentBase, byte baseQuality) {
        if(observedBase == fragmentBase)
            return QualityUtils.MINUS_QUALITY_PROB_LOG10[baseQuality];
        else
            return QualityUtils.QUALITY_PROB_LOG10[baseQuality] - log10_3;
    }

    protected int getFilteredDepth(Pileup pileup) {
        int count = 0;
        for ( PileupReadInfo p : pileup.getTotalPileup() ) {
            if ( !p.isDeletionBase() && BaseUtils.isRegularBase( p.getByteBase() ) )
                count += 1;
        }

        return count;
    }

}

