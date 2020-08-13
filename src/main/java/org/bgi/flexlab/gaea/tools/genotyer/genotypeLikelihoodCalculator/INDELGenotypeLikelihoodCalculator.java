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
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.mapreduce.genotyper.GenotyperOptions;
import org.bgi.flexlab.gaea.util.BaseUtils;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by zhangyong on 2016/12/29.
 *
 * mainly came from GATK 2.3.9-lite
 */
public class INDELGenotypeLikelihoodCalculator extends GenotypeLikelihoodCalculator {
    private static final int HAPLOTYPE_SIZE = 80;

    public static int REF_WIN_Extend = 200;


    private boolean ignoreSNPAllelesWhenGenotypingIndels = false;

    private PairHMMIndelErrorModel pairModel;


    private LinkedHashMap<Allele, Haplotype> haplotypeMap;

    private List<Allele> alleleList = new ArrayList<>();

    /**
     * constrcutor
     * @param options
     */
    public INDELGenotypeLikelihoodCalculator(GenotyperOptions options) {
        super(options);
        haplotypeMap = new LinkedHashMap<>();
        pairModel = new PairHMMIndelErrorModel(options.getIndelGapOpenPenalty(), options.getIndelGapContinuationPenalty(), options.getPairHmmImplementation());
        ignoreSNPAllelesWhenGenotypingIndels = false;
    }

    /**
     * calculate genotype likelihoods for SNP
     *
     * @param mpileup   multi-sample pileups
     * @param reference reference
     * @param options   options
     * @return variant context with likelihoods
     */
    public VariantContext genotypeLikelihoodCalculate(Mpileup mpileup, ChromosomeInformationShare reference, GenotyperOptions options, GenomeLocationParser locationParser, Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        Map<String, Pileup> pileups = mpileup.getCurrentPosPileup();
        int position = mpileup.getPosition();
        if (pileups == null || pileups.isEmpty())
            return null;

        haplotypeMap.clear();
        perReadAlleleLikelihoodMap.clear();

        int winStart = Math.max( position - REF_WIN_Extend, 0 );
        int winStop = Math.min( position + REF_WIN_Extend, reference.getLength() - 1 );
        GenomeLocation refWindows = locationParser.createGenomeLocation(reference.getChromosomeName(), winStart, winStop);

        //construct haplotypes
        ///get alleles from pileup
        //System.err.println("get indel alleles");
        alleleList = getConsensusAlleles(pileups, position, reference, options, locationParser);
        if (alleleList.isEmpty())
            return null;
        ///construct haplotypes
        //System.err.println("construct haplotypes");
        getHaplotypeMapFromAlleles(alleleList, reference, position, refWindows,locationParser, haplotypeMap);
        if (haplotypeMap == null || haplotypeMap.isEmpty())
            return null;

        //System.err.println("make variant context");
        // start making the VariantContext
        // For all non-snp VC types, VC end location is just startLocation + length of ref allele including padding base.
        final int endLoc = position + alleleList.get(0).length() - 1;
        final int eventLength = getEventLength(alleleList);

        final VariantContextBuilder builder = new VariantContextBuilder("UG_call", reference.getChromosomeName(), position + 1, endLoc + 1, alleleList);

        // create the genotypes; no-call everyone for now
        GenotypesContext genotypes = GenotypesContext.create();
        final List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        // For each sample, get genotype likelihoods based on pileup
        // compute prior likelihoods on haplotypes, and initialize haplotype likelihood matrix with them.
        //System.err.println("cal likelihoods");
        for (String sample : pileups.keySet()) {
            Pileup pileup = pileups.get(sample);

            if (!perReadAlleleLikelihoodMap.containsKey(sample)){
                // no likelihoods have been computed for this sample at this site
                PerReadAlleleLikelihoodMap pr = PerReadAlleleLikelihoodMap.getBestAvailablePerReadAlleleLikelihoodMap();
                pr.setRefAllele(alleleList.get(0));
                perReadAlleleLikelihoodMap.put(sample, pr);
            }
            if (pileup != null) {
                final GenotypeBuilder b = new GenotypeBuilder(sample);
                final double[] genotypeLikelihoods = pairModel.computeDiploidReadHaplotypeLikelihoods(pileup, haplotypeMap,
                        reference, eventLength, perReadAlleleLikelihoodMap.get(sample), refWindows, position, options.getContaminationFraction());
                b.PL(genotypeLikelihoods);
                b.DP(getFilteredDepth(pileup));
                genotypes.add(b.make());

            }
        }
        return builder.genotypes(genotypes).make();
    }

    public static List<Allele> getConsensusAlleles(Map<String, Pileup> pileups, int positon, ChromosomeInformationShare reference, GenotyperOptions options, GenomeLocationParser locationParser) {
        ConsensusAlleleCounter counter = new ConsensusAlleleCounter(true, options.getMinIndelCountForGenotyping(), options.getMinIndelFractionPerSample());
        return counter.computeConsensusAlleles(reference, pileups, positon, locationParser);
    }

    public static void getHaplotypeMapFromAlleles(final List<Allele> alleleList,
                                                  final ChromosomeInformationShare ref,
                                                  final int position,
                                                  final GenomeLocation refWindows,
                                                  final GenomeLocationParser locationParser,
                                                  final LinkedHashMap<Allele, Haplotype> haplotypeMap) {
        // protect against having an indel too close to the edge of a contig
        if (position <= HAPLOTYPE_SIZE) {
            haplotypeMap.clear();
        }
        // check if there is enough reference window to create haplotypes (can be an issue at end of contigs)
        else if (ref.getLength() - 1 < position + HAPLOTYPE_SIZE) {
            haplotypeMap.clear();
        }
        else if (alleleList.isEmpty())
            haplotypeMap.clear();
        else {
            final int eventLength = getEventLength(alleleList);
            final int hsize = refWindows.size() - Math.abs(eventLength) - 1;
            final int numPrefBases = position - refWindows.getStart() + 1;
            //System.err.println("event length:" + eventLength + "\thsize:" + hsize + "\tnumPrefBases:" + numPrefBases);
            if (hsize <= 0) { // protect against event lengths larger than ref window sizes
                haplotypeMap.clear();
            }
            else {
                haplotypeMap.putAll(Haplotype.makeHaplotypeListFromAlleles(alleleList, position, ref, refWindows, locationParser, hsize, numPrefBases));
            }
        }
    }

    public static int getEventLength(List<Allele> alleleList) {
        Allele refAllele = alleleList.get(0);
        Allele altAllele = alleleList.get(1);
        // look for alt allele that has biggest length distance to ref allele
        int maxLenDiff = 0;
        for (Allele a : alleleList) {
            if (a.isNonReference()) {
                int lenDiff = Math.abs(a.getBaseString().length() - refAllele.getBaseString().length());
                if (lenDiff > maxLenDiff) {
                    maxLenDiff = lenDiff;
                    altAllele = a;
                }
            }
        }

        return altAllele.getBaseString().length() - refAllele.getBaseString().length();
    }

    // Overload function in GenotypeLikelihoodsCalculationModel so that, for an indel case, we consider a deletion as part of the pileup,
    // so that per-sample DP will include deletions covering the event.
    protected int getFilteredDepth(Pileup pileup) {
        int count = 0;
        for (PileupReadInfo p : pileup.getTotalPileup()) {
            if (p.isDeletionBase() || p.isInsertionAtBeginningOfRead() || BaseUtils.isRegularBase(p.getByteBase()))
                count += 1; //p.getRepresentativeCount();
        }

        return count;
    }
}
