package org.bgi.flexlab.gaea.tools.markduplicate;

import htsjdk.samtools.SAMRecord;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.bgi.flexlab.gaea.util.RandomUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;


public class MarkDuplicatesFunc {
    public void markDup(ArrayList<SAMRecord> sams) {
        Map<String, ReadEnds> readsEnds = new HashMap<>();
        Map<String, ArrayList<SAMRecord>> samPairs = new HashMap<>();

        for (SAMRecord sam : sams){
            ReadEnds ends;
            ArrayList<SAMRecord> samPair;
            if(!readsEnds.containsKey(sam.getReadName())) {
                ends = new ReadEnds();
                buildReadEnds(sam, ends);
                readsEnds.put(sam.getReadName(), ends);

                samPair = new ArrayList<>();
                samPair.add(sam);
                samPairs.put(sam.getReadName(), samPair);
            } else {
                ends = readsEnds.get(sam.getReadName());
                buildReadEnds(sam, ends);

                samPair = samPairs.get(sam.getReadName());
                samPair.add(sam);
            }

        }

        //get clusters
        Map<ClusterIndex, ArrayList<String>> clusters = new HashMap<>();
        for(String rName : readsEnds.keySet()) {
            ReadEnds ends = readsEnds.get(rName);
            ClusterIndex index = new ClusterIndex(ends);
            if(clusters.containsKey(index)) {
                clusters.get(index).add(rName);
            } else {
                ArrayList<String> pairName = new ArrayList<>();
                pairName.add(rName);
                clusters.put(index, pairName);
            }
        }

        // deal cluster
        for (ArrayList<String> rNames: clusters.values()) {

            if(rNames.size() == 1)
                continue;

            //get max score pair
            String maxScoreName = "";
            int maxScore = 0;
            for(String rName : rNames) {
                ReadEnds ends = readsEnds.get(rName);
                //System.err.println("rName:" + rName + "\tscore:" + ends.score);
                if(ends.score >= maxScore) {
                    //System.err.println("old maxScore:" + maxScore + "\tnew maxScore:" + ends.score);
                    //System.err.println("old maxScoreName:" + maxScoreName + "\tnew maxScoreName:" + rName);
                    maxScore = ends.score;
                    maxScoreName = rName;
                }
            }

            //mark rest pairs
            for(String rName : rNames) {
                ArrayList<SAMRecord> samPair = samPairs.get(rName);
                assert !maxScoreName.equals("");
                if(!rName.equals(maxScoreName)){
                    for(SAMRecord sam : samPair) {
                        sam.setDuplicateReadFlag(true);
                    }
                }
            }

        }
    }

    /**
     * build readEnds
     * @param rec
     * @param ends
     */
    private void buildReadEnds(SAMRecord rec, ReadEnds ends) {
        if(ends.read1SequenceIndex == -1) {
            ends.read1SequenceIndex = rec.getReferenceIndex();
            ends.read1Coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();;
            ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
            ends.MQ = rec.getMappingQuality();
            ends.score = getScore(rec);
            // Doing this lets the ends object know that it's part of a pair
            if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
                ends.read2SequenceIndex = rec.getMateReferenceIndex();
            }
        } else {
            //PE
            int sequence = rec.getReferenceIndex();
            int coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
            if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
                if(sequence > ends.read1SequenceIndex || (sequence == ends.read1SequenceIndex && coordinate >= ends.read1Coordinate)) {
                    ends.read2SequenceIndex=sequence;
                    ends.read2Coordinate  = coordinate;
                    ends.orientation = getOrientationByte(ends.orientation == ReadEnds.R, rec.getReadNegativeStrandFlag());
                } else {
                    ends.read2SequenceIndex=ends.read1SequenceIndex;
                    ends.read2Coordinate  = ends.read1Coordinate;
                    ends.read2Index=ends.read1Index;
                    ends.read1SequenceIndex=sequence;
                    ends.read1Coordinate  = coordinate;
                    ends.orientation = getOrientationByte(rec.getReadNegativeStrandFlag(), ends.orientation == ReadEnds.R);
                }
                ends.score += getScore(rec);
                ends.MQ +=rec.getMappingQuality();
            }
        }
    }

    /**
     * get score of a read, maybe we could add MQ into account
     * @param rec
     * @return
     */
    private short getScore(final SAMRecord rec) {
        short score = 0;
        for (final byte b : rec.getBaseQualities()) {
            if (b >= 15) score += b;
        }
        return score;
    }

    /**
     * get pair orientation
     * @param read1NegativeStrand
     * @param read2NegativeStrand
     * @return
     */
    private byte getOrientationByte(final boolean read1NegativeStrand, final boolean read2NegativeStrand) {
        if (read1NegativeStrand) {
            if (read2NegativeStrand)  return ReadEnds.RR;
            else return ReadEnds.RF;
        }
        else {
            if (read2NegativeStrand)  return ReadEnds.FR;
            else return ReadEnds.FF;
        }
    }

    /**
     * class for cluster key
     * @author zy1905
     *
     */
    private class ClusterIndex implements Comparable<ClusterIndex>{
        private byte orientation;
        private int read1SequenceIndex=-1;
        private int read1Coordinate   = -1;
        private int read2SequenceIndex=-1;
        private int read2Coordinate   = -1;

        ClusterIndex(ReadEnds ends) {
            orientation = ends.orientation;
            read1SequenceIndex = ends.read1SequenceIndex;
            read1Coordinate = ends.read1Coordinate;
            read2SequenceIndex = ends.read2SequenceIndex;
            read2Coordinate = ends.read2Coordinate;
        }

        @Override
        public int compareTo(ClusterIndex ci) {
            int retval  = read1SequenceIndex - ci.read1SequenceIndex;
            if (retval == 0) retval = read1Coordinate - ci.read1Coordinate;
            if (retval == 0) retval = orientation - ci.orientation;
            if(read2SequenceIndex != -1) {
                if (retval == 0) retval = read2SequenceIndex - ci.read2SequenceIndex;
                if (retval == 0) retval = read2Coordinate - ci.read2Coordinate;
            }
            return retval;
        }

        @Override
        public int hashCode() {
            return new HashCodeBuilder().append(read1SequenceIndex).append(read1Coordinate).append(orientation).toHashCode();
        }

        @Override
        public boolean equals(Object obj) {
            if(!(obj instanceof ClusterIndex)) {
                return false;
            } else if(obj == this)
                return true;

            ClusterIndex index = (ClusterIndex) obj;
            return this.compareTo(index) == 0;
        }
    }

}
