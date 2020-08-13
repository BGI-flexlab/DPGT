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

package org.bgi.flexlab.gaea.tools.genotyer.annotator;



import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.alignment.AlignmentsBasic;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.pileup.PileupReadInfo;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.StandardAnnotation;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.Haplotype;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;
import org.bgi.flexlab.gaea.util.BaseUtils;
import org.bgi.flexlab.gaea.util.MathUtils;
import org.bgi.flexlab.gaea.util.QualityUtils;

import java.io.Serializable;
import java.util.*;
/**
 * Consistency of the site with two (and only two) segregating haplotypes. Higher scores
 * are indicative of regions with bad alignments, often leading to artifactual SNP and indel calls.
 * Note that the Haplotype Score is only calculated for sites with read coverage.
 */
public class HaplotypeScore extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
    private final static boolean DEBUG = false;
    private final static int MIN_CONTEXT_WING_SIZE = 10;
    private final static int MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER = 50;
    private final static char REGEXP_WILDCARD = '.';

    public Map<String, Object> annotate(final VariantDataTracker tracker,
                                        final ChromosomeInformationShare ref,
                                        final Mpileup mpileup,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
           	if (vc.isSNP() && mpileup != null)
            return annotatePileup(ref, mpileup, vc);
        else if (stratifiedPerReadAlleleLikelihoodMap != null && vc.isVariant())
            return annotateWithLikelihoods(stratifiedPerReadAlleleLikelihoodMap, vc);
        else
            return null;
    }

    private Map<String, Object> annotatePileup(final ChromosomeInformationShare ref,
                                        final Mpileup mpileup,
                                        final VariantContext vc) {

        if (mpileup.getSize() == 0) // size 0 means that call was made by someone else and we have no data here
            return null;

        final int contextWingSize = Math.min((400 - 1) / 2, MIN_CONTEXT_WING_SIZE);
        final int contextSize = contextWingSize * 2 + 1;

        final int locus = mpileup.getPosition();

        // Compute all haplotypes consistent with the current read pileup
        final List<Haplotype> haplotypes = computeHaplotypes(mpileup.joinPileups(), contextSize, locus, vc);

        final MathUtils.RunningAverage scoreRA = new MathUtils.RunningAverage();
        if (haplotypes != null) {
            for (final Genotype genotype : vc.getGenotypes()) {
                final Pileup thisPileup = mpileup.getCurrentPosPileup().get(genotype.getSampleName());
                if (thisPileup != null) {
                    scoreRA.add(scoreReadsAgainstHaplotypes(haplotypes, thisPileup, contextSize, locus)); // Taking the simple average of all sample's score since the score can be negative and the RMS doesn't make sense
                }
            }
        }

        // annotate the score in the info field
        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.4f", scoreRA.mean()));
        return map;
    }

    private Map<String, Object> annotateWithLikelihoods(final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                        final VariantContext vc) {


        final MathUtils.RunningAverage scoreRA = new MathUtils.RunningAverage();
        for (final Genotype genotype : vc.getGenotypes()) {
            final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());
            if (perReadAlleleLikelihoodMap == null)
                continue;

            Double d = scoreIndelsAgainstHaplotypes(perReadAlleleLikelihoodMap);
            if (d == null)
                continue;
           scoreRA.add(d); // Taking the simple average of all sample's score since the score can be negative and the RMS doesn't make sense
        }

 //       if (scoreRA.observationCount() == 0)
 //           return null;

        // annotate the score in the info field
        final Map<String, Object> map = new HashMap<String, Object>();
        map.put(getKeyNames().get(0), String.format("%.4f", scoreRA.mean()));
        return map;

    }

    private static class HaplotypeComparator implements Comparator<Haplotype>, Serializable {

        public int compare(Haplotype a, Haplotype b) {
            if (a.getQualitySum() < b.getQualitySum())
                return 1;
            if (a.getQualitySum() > b.getQualitySum()) {
                return -1;
            }
            return 0;
        }
    }

    private List<Haplotype> computeHaplotypes(final Pileup pileup, final int contextSize, final int locus, final VariantContext vc) {
        // Compute all possible haplotypes consistent with current pileup

        int haplotypesToCompute = vc.getAlternateAlleles().size() + 1;

        final PriorityQueue<Haplotype> candidateHaplotypeQueue = new PriorityQueue<Haplotype>(100, new HaplotypeComparator());
        final PriorityQueue<Haplotype> consensusHaplotypeQueue = new PriorityQueue<Haplotype>(MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER, new HaplotypeComparator());

        for (final PileupReadInfo p : pileup.getTotalPileup()) {
            final Haplotype haplotypeFromRead = getHaplotypeFromRead(p, contextSize, locus);
            candidateHaplotypeQueue.add(haplotypeFromRead);
        }

        // Now that priority queue has been built with all reads at context, we need to merge and find possible segregating haplotypes
        Haplotype elem;
        while ((elem = candidateHaplotypeQueue.poll()) != null) {
            boolean foundHaplotypeMatch = false;
            Haplotype lastCheckedHaplotype = null;
            for (final Haplotype haplotypeFromList : consensusHaplotypeQueue) {
                final Haplotype consensusHaplotype = getConsensusHaplotype(elem, haplotypeFromList);
                if (consensusHaplotype != null) {
                    foundHaplotypeMatch = true;
                    if (consensusHaplotype.getQualitySum() > haplotypeFromList.getQualitySum()) {
                        consensusHaplotypeQueue.remove(haplotypeFromList);
                        consensusHaplotypeQueue.add(consensusHaplotype);
                    }
                    break;
                } else {
                    lastCheckedHaplotype = haplotypeFromList;
                }
            }

            if (!foundHaplotypeMatch && consensusHaplotypeQueue.size() < MAX_CONSENSUS_HAPLOTYPES_TO_CONSIDER) {
                consensusHaplotypeQueue.add(elem);
            } else if (!foundHaplotypeMatch && lastCheckedHaplotype != null && elem.getQualitySum() > lastCheckedHaplotype.getQualitySum()) {
                consensusHaplotypeQueue.remove(lastCheckedHaplotype);
                consensusHaplotypeQueue.add(elem);
            }
        }

        // Now retrieve the N most popular haplotypes
        if (consensusHaplotypeQueue.size() > 0) {
            // The consensus haplotypes are in a quality-ordered priority queue, so the best haplotypes are just the ones at the front of the queue
            final Haplotype haplotype1 = consensusHaplotypeQueue.poll();

            List<Haplotype> hlist = new ArrayList<Haplotype>();
            hlist.add(new Haplotype(haplotype1.getBases(), 60));

            for (int k = 1; k < haplotypesToCompute; k++) {
                Haplotype haplotype2 = consensusHaplotypeQueue.poll();
                if (haplotype2 == null) {
                    haplotype2 = haplotype1;
                } // Sometimes only the reference haplotype can be found
                hlist.add(new Haplotype(haplotype2.getBases(), 20));
            }
            return hlist;
        } else
            return null;
    }

    private Haplotype getHaplotypeFromRead(final PileupReadInfo p, final int contextSize, final int locus) {
        final AlignmentsBasic read = p.getReadInfo();

        final byte[] haplotypeBases = new byte[contextSize];
        Arrays.fill(haplotypeBases, (byte) REGEXP_WILDCARD);
        final double[] baseQualities = new double[contextSize];
        Arrays.fill(baseQualities, 0.0);

        byte[] readBases = read.getReadBases();
        readBases = p.readToAlignmentByteArray(read.getCigars(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getQualities();
        readQuals = p.readToAlignmentByteArray(read.getCigars(), readQuals); // Shift the location of the qual scores based on the Cigar string

        final int readOffsetFromPileup = p.calcAlignmentByteArrayOffset(read.getPosition(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1) / 2;

        for (int i = 0; i < contextSize; i++) {
            final int baseOffset = i + baseOffsetStart;
            if (baseOffset < 0) {
                continue;
            }
            if (baseOffset >= readBases.length) {
                break;
            }
            if (readQuals[baseOffset] == PileupReadInfo.DELETION_BASE) {
                readQuals[baseOffset] = PileupReadInfo.DELETION_QUAL;
            }
            if (!BaseUtils.isRegularBase(readBases[baseOffset])) {
                readBases[baseOffset] = (byte) REGEXP_WILDCARD;
                readQuals[baseOffset] = (byte) 0;
            } // N's shouldn't be treated as distinct bases
            readQuals[baseOffset] = (byte) Math.min((int) readQuals[baseOffset], p.getMappingQuality());
            if (((int) readQuals[baseOffset]) < 5) {
                readQuals[baseOffset] = (byte) 0;
            } // quals less than 5 are used as codes and don't have actual probabilistic meaning behind them
            haplotypeBases[i] = readBases[baseOffset];
            baseQualities[i] = (double) readQuals[baseOffset];
        }

        return new Haplotype(haplotypeBases, baseQualities);
    }

    private Haplotype getConsensusHaplotype(final Haplotype haplotypeA, final Haplotype haplotypeB) {
        final byte[] a = haplotypeA.getBases();
        final byte[] b = haplotypeB.getBases();

        if (a.length != b.length) {
            throw new UserException("Haplotypes a and b must be of same length");
        }

        byte chA, chB;
        final byte wc = (byte) REGEXP_WILDCARD;

        final int length = a.length;
        final byte[] consensusChars = new byte[length];
        final double[] consensusQuals = new double[length];

        final double[] qualsA = haplotypeA.getQuals();
        final double[] qualsB = haplotypeB.getQuals();

        for (int i = 0; i < length; i++) {
            chA = a[i];
            chB = b[i];

            if ((chA != chB) && (chA != wc) && (chB != wc))
                return null;

            if ((chA == wc) && (chB == wc)) {
                consensusChars[i] = wc;
                consensusQuals[i] = 0.0;
            } else if ((chA == wc)) {
                consensusChars[i] = chB;
                consensusQuals[i] = qualsB[i];
            } else if ((chB == wc)) {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i];
            } else {
                consensusChars[i] = chA;
                consensusQuals[i] = qualsA[i] + qualsB[i];
            }
        }

        return new Haplotype(consensusChars, consensusQuals);
    }

    // calculate the haplotype scores by walking over all reads and comparing them to the haplotypes
    private double scoreReadsAgainstHaplotypes(final List<Haplotype> haplotypes, final Pileup pileup, final int contextSize, final int locus) {
        if (DEBUG) System.out.printf("HAP1: %s%n", haplotypes.get(0));
        if (DEBUG) System.out.printf("HAP2: %s%n", haplotypes.get(1));

        final ArrayList<double[]> haplotypeScores = new ArrayList<double[]>();
        for (final PileupReadInfo p : pileup.getTotalPileup()) {
            // Score all the reads in the pileup, even the filtered ones
            final double[] scores = new double[haplotypes.size()];
            for (int i = 0; i < haplotypes.size(); i++) {
                final Haplotype haplotype = haplotypes.get(i);
                final double score = scoreReadAgainstHaplotype(p, contextSize, haplotype, locus);
                scores[i] = score;                
            }
            haplotypeScores.add(scores);
        }

        double overallScore = 0.0;
        for (final double[] readHaplotypeScores : haplotypeScores) {
            overallScore += MathUtils.arrayMin(readHaplotypeScores);
        }

        return overallScore;
    }

    private double scoreReadAgainstHaplotype(final PileupReadInfo p, final int contextSize, final Haplotype haplotype, final int locus) {
        double expected = 0.0;
        double mismatches = 0.0;

        // What's the expected mismatch rate under the model that this read is actually sampled from
        // this haplotype?  Let's assume the consensus base c is a random choice one of A, C, G, or T, and that
        // the observed base is actually from a c with an error rate e.  Since e is the rate at which we'd
        // see a miscalled c, the expected mismatch rate is really e.  So the expected number of mismatches
        // is just sum_i e_i for i from 1..n for n sites
        //
        // Now, what's the probabilistic sum of mismatches?  Suppose that the base b is equal to c.  Well, it could
        // actually be a miscall in a matching direction, which would happen at a e / 3 rate.  If b != c, then
        // the chance that it is actually a mismatch is 1 - e, since any of the other 3 options would be a mismatch.
        // so the probability-weighted mismatch rate is sum_i ( matched ? e_i / 3 : 1 - e_i ) for i = 1 ... n
        final byte[] haplotypeBases = haplotype.getBases();
        final AlignmentsBasic read = p.getReadInfo();
        byte[] readBases = read.getReadBases();

        readBases = p.readToAlignmentByteArray(p.getReadInfo().getCigars(), readBases); // Adjust the read bases based on the Cigar string
        byte[] readQuals = read.getQualities();
        readQuals = p.readToAlignmentByteArray(p.getReadInfo().getCigars(), readQuals); // Shift the location of the qual scores based on the Cigar string
        int readOffsetFromPileup = p.getQpos();
        readOffsetFromPileup = p.calcAlignmentByteArrayOffset(read.getPosition(), locus);
        final int baseOffsetStart = readOffsetFromPileup - (contextSize - 1) / 2;

        for (int i = 0; i < contextSize; i++) {
            final int baseOffset = i + baseOffsetStart;
            if (baseOffset < 0) {
                continue;
            }
            if (baseOffset >= readBases.length) {
                break;
            }

            final byte haplotypeBase = haplotypeBases[i];
            final byte readBase = readBases[baseOffset];

            final boolean matched = (readBase == haplotypeBase || haplotypeBase == (byte) REGEXP_WILDCARD);
            byte qual = readQuals[baseOffset];
            if (qual == PileupReadInfo.DELETION_BASE) {
                qual = PileupReadInfo.DELETION_QUAL;
            } // calcAlignmentByteArrayOffset fills the readQuals array with DELETION_BASE at deletions
            qual = (byte) Math.min((int) qual, p.getMappingQuality());
            if (((int) qual) >= 5) { // quals less than 5 are used as codes and don't have actual probabilistic meaning behind them
                final double e = QualityUtils.qualityToErrorProbability(qual);
                expected += e;
                mismatches += matched ? e : 1.0 - e / 3.0;
            }

            // a more sophisticated calculation would include the reference quality, but it's nice to actually penalize
            // the mismatching of poorly determined regions of the consensus
        }

        return mismatches - expected;
    }


    private Double scoreIndelsAgainstHaplotypes(final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap) {
        final ArrayList<double[]> haplotypeScores = new ArrayList<double[]>();

        if (perReadAlleleLikelihoodMap.isEmpty())
            return null;

        for (Map<Allele,Double> el : perReadAlleleLikelihoodMap.getLikelihoodMapValues()) {

            // retrieve likelihood information corresponding to this read
            // Score all the reads in the pileup, even the filtered ones
            final double[] scores = new double[el.size()];
            int i = 0;
            for (Map.Entry<Allele, Double> a : el.entrySet()) {
                scores[i++] = -a.getValue();                
            }

            haplotypeScores.add(scores);
        }

        // indel likelihoods are strict log-probs, not phred scored
        double overallScore = 0.0;
        for (final double[] readHaplotypeScores : haplotypeScores) {
            overallScore += MathUtils.arrayMin(readHaplotypeScores);
        }

        return overallScore;

    }


    public List<String> getKeyNames() {
        return Arrays.asList("HaplotypeScore");
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFInfoHeaderLine("HaplotypeScore", 1, VCFHeaderLineType.Float, "Consistency of the site with at most two segregating haplotypes"));
    }
}
