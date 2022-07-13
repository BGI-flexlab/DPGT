package org.bgi.flexlab.dpgt.jointcalling;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class VariantSiteSet {
    private static final Logger logger = LoggerFactory.getLogger(VariantSiteSet.class);
    public BitSet data;
    public int start;
    public int end;
    public int size;

    public VariantSiteSet(int start, int end) {
        this.start = start;
        this.end = end;
        this.size = end - start + 1;
        this.data = new BitSet(this.size);
    }

    public VariantSiteSet(int start, int end, BitSet data) {
        this.start = start;
        this.end = end;
        this.size = end - start + 1;
        this.data = data;
    }

    void set(int pos) {
        validatePosition(pos);
        data.set(pos - this.start);
    }

    void reset(int pos) {
        validatePosition(pos);
        data.set(pos - this.start, false);
    }

    boolean test(int pos) {
        validatePosition(pos);
        return data.get(pos - this.start);
    }

    boolean test(int start, int end) {
        validatePosition(start);
        validatePosition(end);
        return !data.get(start - this.start, end - this.start).isEmpty();
    }

    boolean get(int pos) {
        validatePosition(pos);
        return data.get(pos - this.start);
    }

    BitSet get(int start, int end) {
        validatePosition(start);
        validatePosition(end);
        return data.get(start - this.start, end - this.start);
    }

    private void validatePosition(int pos) {
        if (pos < this.start || pos > this.end) {
            logger.error("Position %l out of range [%l, %l]", pos, start, end);
            System.exit(1);
        }
    }

}
