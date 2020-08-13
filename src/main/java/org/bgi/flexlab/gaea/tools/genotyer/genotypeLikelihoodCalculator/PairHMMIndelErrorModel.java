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

import htsjdk.variant.variantcontext.Allele;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.bam.clipper.ReadClipper;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.util.MathUtils;
import org.bgi.flexlab.gaea.util.ReadUtils;
import org.bgi.flexlab.gaea.util.pairhmm.ExactPairHMM;
import org.bgi.flexlab.gaea.util.pairhmm.OriginalPairHMM;
import org.bgi.flexlab.gaea.util.pairhmm.PairHMM;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * mainly came from GATK 2.3.9-lite
 */
public class PairHMMIndelErrorModel {
    public static final int BASE_QUAL_THRESHOLD = 20;

    private static final int MAX_CACHED_QUAL = 127;

    private static final double baseMatchArray[];
    private static final double baseMismatchArray[];

    private final static double LOG_ONE_HALF;

    private static final int START_HRUN_GAP_IDX = 4;
    private static final int MAX_HRUN_GAP_IDX = 20;

    private static final byte MIN_GAP_OPEN_PENALTY = 30;
    private static final byte MIN_GAP_CONT_PENALTY = 10;
    private static final byte GAP_PENALTY_HRUN_STEP = 1; // each increase in hrun decreases gap penalty by this.

    private final byte[] GAP_OPEN_PROB_TABLE;
    private final byte[] GAP_CONT_PROB_TABLE;

    private final PairHMM pairHMM;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////

    static {
        LOG_ONE_HALF= -Math.log10(2.0);

        baseMatchArray = new double[MAX_CACHED_QUAL+1];
        baseMismatchArray = new double[MAX_CACHED_QUAL+1];
        for (int k=1; k <= MAX_CACHED_QUAL; k++) {
            double baseProb = Math.pow(10, -k/10.);


            baseMatchArray[k] =  Math.log10(1-baseProb);
            baseMismatchArray[k] = Math.log10(baseProb);
        }
    }

    public PairHMMIndelErrorModel(byte indelGOP, byte indelGCP, final PairHMM.HMM_IMPLEMENTATION hmmType ) {

        switch (hmmType) {
            case EXACT:
                pairHMM = new ExactPairHMM();
                break;
            case ORIGINAL:
                pairHMM = new OriginalPairHMM();
                break;
            case CACHING:
            case LOGLESS_CACHING:
            default:
                throw new UserException.BadArgumentValueException("pairHMM", "Specified pairHMM implementation is unrecognized or incompatible with the UnifiedGenotyper. Acceptable options are ORIGINAL and EXACT.");
        }

        // fill gap penalty table, affine naive model:
        this.GAP_CONT_PROB_TABLE = new byte[MAX_HRUN_GAP_IDX];
        this.GAP_OPEN_PROB_TABLE = new byte[MAX_HRUN_GAP_IDX];

        for (int i = 0; i < START_HRUN_GAP_IDX; i++) {
            GAP_OPEN_PROB_TABLE[i] = indelGOP;
            GAP_CONT_PROB_TABLE[i] = indelGCP;
        }

        double step = GAP_PENALTY_HRUN_STEP/10.0;

        // initialize gop and gcp to their default values
        byte gop = indelGOP;
        byte gcp = indelGCP;

        // all of the following is computed in QUal-space
        for (int i=START_HRUN_GAP_IDX; i < MAX_HRUN_GAP_IDX; i++) {
            gop -= GAP_PENALTY_HRUN_STEP;
            if (gop < MIN_GAP_OPEN_PENALTY)
                gop = MIN_GAP_OPEN_PENALTY;

            gcp -= step;
            if(gcp < MIN_GAP_CONT_PENALTY)
                gcp = MIN_GAP_CONT_PENALTY;
            GAP_OPEN_PROB_TABLE[i] = gop;
            GAP_CONT_PROB_TABLE[i] = gcp;
        }

    }

    static private void getContextHomopolymerLength(final byte[] refBytes, final int[] hrunArray) {
        // compute forward hrun length, example:
        // AGGTGACCCCCCTGAGAG
        // 001000012345000000
        hrunArray[0] = 0;
        int[] hforward = new int[hrunArray.length];
        int[] hreverse = new int[hrunArray.length];

        for (int i = 1; i < refBytes.length; i++) {
            if (refBytes[i] == refBytes[i-1])
                hforward[i] = hforward[i-1]+1;
            else
                hforward[i] = 0;
        }

        // do similar thing for reverse length, example:
        // AGGTGACCCCCCTGAGAG
        // 021000543210000000
        // and then accumulate with forward values.
        // Total:
        // AGGTGACCCCCCTGAGAG
        // 022000555555000000
        for (int i=refBytes.length-1; i > 0; i--) {
            if (refBytes[i-1] == refBytes[i])
                hreverse[i-1] += hreverse[i]+1;
        }

        for (int i = 1; i < refBytes.length; i++)
            hrunArray[i] = hforward[i]+hreverse[i];
    }


    private void fillGapProbabilities(final int[] hrunProfile,
                                      final byte[] contextLogGapOpenProbabilities,
                                      final byte[] contextLogGapContinuationProbabilities) {
        // fill based on lookup table
        for (int i = 0; i < hrunProfile.length; i++) {
            if (hrunProfile[i] >= MAX_HRUN_GAP_IDX) {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[MAX_HRUN_GAP_IDX-1];
            }
            else {
                contextLogGapOpenProbabilities[i] = GAP_OPEN_PROB_TABLE[hrunProfile[i]];
                contextLogGapContinuationProbabilities[i] = GAP_CONT_PROB_TABLE[hrunProfile[i]];
            }
        }
    }


    public synchronized double[] computeDiploidReadHaplotypeLikelihoods(final Pileup pileup,
                                                                        final LinkedHashMap<Allele, Haplotype> haplotypeMap,
                                                                        final ChromosomeInformationShare ref,
                                                                        final int eventLength,
                                                                        final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap,
                                                                        final GenomeLocation refWindows,
                                                                        final int position,
                                                                        final double downsamplingFraction
                                                                        ) {
        final int numHaplotypes = haplotypeMap.size();

        final int readCounts[] = new int[pileup.getNumberOfElements()];
        final double[][] readLikelihoods = computeGeneralReadHaplotypeLikelihoods(pileup, haplotypeMap, ref, eventLength, perReadAlleleLikelihoodMap, refWindows, position, readCounts);
        perReadAlleleLikelihoodMap.performPerAlleleDownsampling(downsamplingFraction);
        return getDiploidHaplotypeLikelihoods(numHaplotypes, readCounts, readLikelihoods);
        
    }

    //@Ensures("result != null && result.length == pileup.getNumberOfElements()")
    public synchronized double[][] computeGeneralReadHaplotypeLikelihoods(final Pileup pileup,
                                                                          final LinkedHashMap<Allele, Haplotype> haplotypeMap,
                                                                          final ChromosomeInformationShare ref,
                                                                          final int eventLength, 
                                                                          final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap,
                                                                          final GenomeLocation refWindows,
                                                                          final int position,
                                                                          final int[] readCounts) {
        final double readLikelihoods[][] = new double[pileup.getNumberOfElements()][haplotypeMap.size()];

        int readIdx=0;
        for (PileupReadInfo p: pileup.getFilteredPileup()) {
            // > 1 when the read is a consensus read representing multiple independent observations
            readCounts[readIdx] = 1; //p.getRepresentativeCount();

            // check if we've already computed likelihoods for this pileup element (i.e. for this read at this location)
            if (perReadAlleleLikelihoodMap.containsPileupElement(p)) {
                Map<Allele,Double> el = perReadAlleleLikelihoodMap.getLikelihoodsAssociatedWithPileupElement(p);
                int j=0;
                for (Allele a: haplotypeMap.keySet()) {
                    readLikelihoods[readIdx][j++] = el.get(a);
                }
            }
            else {
                final int refWindowStart = refWindows.getStart();
                final int refWindowStop  = refWindows.getStop();

                //FIXME:: maybe need to be fixed
                //GaeaSAMRecord read = ReadClipper.hardClipAdaptorSequence(p.getRead());
                AlignmentsBasic read = p.getReadInfo();

                read.checkCigar();
                int readSoftStart = read.getSoftStart();
                int readSoftEnd = read.getSoftEnd();
                //System.err.println("read soft start:" + readSoftStart + "\tread soft end:" + readSoftEnd);
                if (!read.isEmpty() && (readSoftEnd > refWindowStop && readSoftStart < refWindowStop))
                    read = ReadClipper.hardClipByReferenceCoordinatesRightTail(read, refWindows.getStop());

                if (!read.isEmpty() && (readSoftStart < refWindowStart && readSoftEnd > refWindowStart))
                    read = ReadClipper.hardClipByReferenceCoordinatesLeftTail (read, refWindows.getStart());

                if (read.isEmpty())
                    continue;

                read.checkCigar();
                // hard-clip low quality ends - this may introduce extra H elements in CIGAR string
                read = ReadClipper.hardClipLowQualEnds(read, (byte) BASE_QUAL_THRESHOLD );

                if (read.isEmpty())
                    continue;

                // get bases of candidate haplotypes that overlap with reads
                final int trailingBases = 3;
                final long readStart = read.getSoftStart();
                final long readEnd = read.getSoftEnd();

                // see if we want to use soft clipped bases. Aligners may soft clip all bases at insertions because they don't match,
                // but they're actually consistent with the insertion!
                // Rule: if a read starts in interval [eventStart-eventLength,eventStart+1] and we are at an insertion, we'll use all soft clipped bases at the beginning.
                // Conversely, if a read ends at [eventStart,eventStart+eventLength] we'll use all soft clipped bases in the end of the read.
                final long eventStartPos = position;

                // compute total number of clipped bases (soft or hard clipped) and only use them if necessary
                final boolean softClips = useSoftClippedBases(read, eventStartPos, eventLength);
                final int numStartSoftClippedBases = softClips ? read.getPosition()- read.getSoftStart() : 0;
                final int numEndSoftClippedBases = softClips ? read.getSoftEnd()- p.getEnd() : 0 ;
                final byte [] unclippedReadBases = read.getReadBases();
                final byte [] unclippedReadQuals = read.getQualities();
                final int extraOffset = Math.abs(eventLength);

                //System.err.println("read:" + read.getPosition());

                /**
                 * Compute genomic locations that candidate haplotypes will span.
                 * Read start and stop locations (variables readStart and readEnd) are the original unclipped positions from SAMRecord,
                 * adjusted by hard clips from Cigar string and by qual-based soft-clipping performed above.
                 * We will propose haplotypes that overlap the read with some padding.
                 * True read start = readStart + numStartSoftClippedBases - ReadUtils.getFirstInsertionOffset(read)
                 * Last term is because if a read starts with an insertion then these bases are not accounted for in readStart.
                 * trailingBases is a padding constant(=3) and we additionally add abs(eventLength) to both sides of read to be able to
                 * differentiate context between two haplotypes
                 */
                long startLocationInRefForHaplotypes = Math.max(readStart + numStartSoftClippedBases - trailingBases - ReadUtils.getFirstInsertionOffset(read)-extraOffset, 0);
                long stopLocationInRefForHaplotypes =  readEnd -numEndSoftClippedBases  + trailingBases + ReadUtils.getLastInsertionOffset(read)+extraOffset;

                int readLength = read.getReadLength()-numStartSoftClippedBases-numEndSoftClippedBases;

                if (startLocationInRefForHaplotypes < refWindows.getStart()) {
                    startLocationInRefForHaplotypes = refWindows.getStart();                                       // read starts before haplotype: read will have to be cut numStartSoftClippedBases += ref.getWindow().getStart() - startLocationInRefForHaplotypes;
                }
                else if (startLocationInRefForHaplotypes > refWindows.getStop()) {
                    startLocationInRefForHaplotypes = refWindows.getStop();                                        // read starts after haplotype: read will have to be clipped completely;
                }

                if (stopLocationInRefForHaplotypes > refWindows.getStop()) {
                    stopLocationInRefForHaplotypes = refWindows.getStop();                                         // check also if end of read will go beyond reference context
                }

                if (stopLocationInRefForHaplotypes <= startLocationInRefForHaplotypes + readLength) {
                    stopLocationInRefForHaplotypes = startLocationInRefForHaplotypes + readLength-1;                    // if there's an insertion in the read, the read stop position will be less than start + read legnth, but we want to compute likelihoods in the whole region that a read might overlap
                }

                // ok, we now figured out the total number of clipped bases on both ends.
                // Figure out where we want to place the haplotype to score read against

                /*System.err.println("softClips:" + softClips + "\tnumStartSoftClippedBases:" + numStartSoftClippedBases +
                        "\tnumEndSoftClippedBases:" + numEndSoftClippedBases + "\textraOffset:" + extraOffset +
                        "\tstartLocationInRefForHaplotypes:" + startLocationInRefForHaplotypes + "\tstopLocationInRefForHaplotypes:"
                        + stopLocationInRefForHaplotypes + "\treadLength:" + readLength);*/

               // LinkedHashMap<Allele,Double> readEl = new LinkedHashMap<Allele,Double>();

                /**
                 * Check if we'll end up with an empty read once all clipping is done
                 */
                if (numStartSoftClippedBases + numEndSoftClippedBases >= unclippedReadBases.length) {
                    int j=0;
                    for (Allele a: haplotypeMap.keySet()) {
                        perReadAlleleLikelihoodMap.add(p,a,0.0);
                        readLikelihoods[readIdx][j++] = 0.0;
                    }
                }
                else {
                    final byte[] readBases = Arrays.copyOfRange(unclippedReadBases,numStartSoftClippedBases, unclippedReadBases.length-numEndSoftClippedBases);
                    final byte[] readQuals = Arrays.copyOfRange(unclippedReadQuals,numStartSoftClippedBases, unclippedReadBases.length-numEndSoftClippedBases);
                    int j=0;

                    byte[] previousHaplotypeSeen = null;
                    final byte[] contextLogGapOpenProbabilities = new byte[readBases.length];
                    final byte[] contextLogGapContinuationProbabilities  = new byte[readBases.length];

                    // get homopolymer length profile for current haplotype
                    final int[] hrunProfile = new int[readBases.length];
                    getContextHomopolymerLength(readBases,hrunProfile);
                    fillGapProbabilities(hrunProfile, contextLogGapOpenProbabilities, contextLogGapContinuationProbabilities);

                    for (Allele a: haplotypeMap.keySet()) {

                        Haplotype haplotype = haplotypeMap.get(a);

                        if (stopLocationInRefForHaplotypes > haplotype.getStopPosition())
                            stopLocationInRefForHaplotypes = haplotype.getStopPosition();

                        if (startLocationInRefForHaplotypes < haplotype.getStartPosition())
                            startLocationInRefForHaplotypes = haplotype.getStartPosition();
                        else if (startLocationInRefForHaplotypes > haplotype.getStopPosition())
                            startLocationInRefForHaplotypes = haplotype.getStopPosition();

                        final long indStart = startLocationInRefForHaplotypes - haplotype.getStartPosition();
                        final long indStop =  stopLocationInRefForHaplotypes - haplotype.getStartPosition();

                        /*System.err.println("startLocationInRefForHaplotypes:" + startLocationInRefForHaplotypes +
                        "\tstopLocationInRefForHaplotypes:" + stopLocationInRefForHaplotypes + "\thaplotype start:" + haplotype.getStartPosition() +
                        "\tindStart:" + indStart + "\tindStop:" + indStop);*/

                        double readLikelihood;

                        final byte[] haplotypeBases = Arrays.copyOfRange(haplotype.getBases(),
                                (int)indStart, (int)indStop);

                        final int X_METRIC_LENGTH = readBases.length+2;
                        final int Y_METRIC_LENGTH = haplotypeBases.length+2;

                        if (previousHaplotypeSeen == null) {
                            //no need to reallocate arrays for each new haplotype, as length won't change
                            pairHMM.initialize(X_METRIC_LENGTH, Y_METRIC_LENGTH);
                        }

                        int startIndexInHaplotype = 0;
                        if (previousHaplotypeSeen != null)
                            startIndexInHaplotype = computeFirstDifferingPosition(haplotypeBases, previousHaplotypeSeen);
                        previousHaplotypeSeen = haplotypeBases.clone();

                        readLikelihood = pairHMM.computeReadLikelihoodGivenHaplotypeLog10(haplotypeBases, readBases, readQuals,
                                //(read.hasBaseIndelQualities() ? read.getBaseInsertionQualities() : contextLogGapOpenProbabilities),
                                //(read.hasBaseIndelQualities() ? read.getBaseDeletionQualities() : contextLogGapOpenProbabilities),
                                contextLogGapOpenProbabilities, contextLogGapOpenProbabilities,
                                contextLogGapContinuationProbabilities, startIndexInHaplotype, false);

                        perReadAlleleLikelihoodMap.add(p, a, readLikelihood);
                        readLikelihoods[readIdx][j++] = readLikelihood;
                    }
                }
            }
            readIdx++;
        }
        return readLikelihoods;
    }

    private boolean useSoftClippedBases(AlignmentsBasic read, long eventStartPos, int eventLength) {
        return !((read.getPosition() >= eventStartPos-eventLength && read.getPosition() <= eventStartPos+1) || (read.getPosition() >= eventStartPos && read.getPosition() <= eventStartPos + eventLength));
    }

    private int computeFirstDifferingPosition(byte[] b1, byte[] b2) {
        if (b1.length != b2.length)
            return 0; // sanity check

        for (int i=0; i < b1.length; i++ ){
            if ( b1[i]!= b2[i] )
                return i;
        }
        return b1.length;
    }

    private static double[] getDiploidHaplotypeLikelihoods(final int numHaplotypes, final int readCounts[], final double readLikelihoods[][]) {
        final double[][] haplotypeLikehoodMatrix = new double[numHaplotypes][numHaplotypes];

        // todo: MAD 09/26/11 -- I'm almost certain this calculation can be simplified to just a single loop without the intermediate NxN matrix
        for (int i=0; i < numHaplotypes; i++) {
            for (int j=i; j < numHaplotypes; j++){
                // combine likelihoods of haplotypeLikelihoods[i], haplotypeLikelihoods[j]
                // L(Hi, Hj) = sum_reads ( Pr(R|Hi)/2 + Pr(R|Hj)/2)
                //readLikelihoods[k][j] has log10(Pr(R_k) | H[j] )
                for (int readIdx = 0; readIdx < readLikelihoods.length; readIdx++) {
                    // Compute log10(10^x1/2 + 10^x2/2) = log10(10^x1+10^x2)-log10(2)
                    // First term is approximated by Jacobian log with table lookup.
                    if (Double.isInfinite(readLikelihoods[readIdx][i]) && Double.isInfinite(readLikelihoods[readIdx][j]))
                        continue;
                    final double li = readLikelihoods[readIdx][i];
                    final double lj = readLikelihoods[readIdx][j];
                    final int readCount = readCounts[readIdx];
                    haplotypeLikehoodMatrix[i][j] += readCount * (MathUtils.approximateLog10SumLog10(li, lj) + LOG_ONE_HALF);
                }
            }
        }

        final double[] genotypeLikelihoods = new double[numHaplotypes*(numHaplotypes+1)/2];
        int k=0;
        for (int j=0; j < numHaplotypes; j++) {
            for (int i=0; i <= j; i++){
                genotypeLikelihoods[k++] = haplotypeLikehoodMatrix[i][j];
            }
        }

        // renormalize so that max element is zero.
        return MathUtils.normalizeFromLog10(genotypeLikelihoods, false, true);
    }
}
