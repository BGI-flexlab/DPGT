package org.bgi.flexlab.gaea.tools.markduplicate;

/**
 * Created by huangzhibo on 2017/4/14.
 */
public class ReadEnds {
    public static final byte F=0, R=1, FF=2, FR=3, RR=4, RF=5;

    public short score = 0;
    public byte orientation;
    public int read1SequenceIndex=-1;
    public int read1Coordinate   = -1;
    public int read1Index=-1;
    public int read2SequenceIndex=-1;
    public int read2Coordinate   = -1;
    public int read2Index=-1;
    public int MQ=-1;
    public int ReadL=-1;
    // Information used to detect optical dupes

    public boolean isPaired() { return !(read2SequenceIndex==-1); }
}
