package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;


import java.io.Serializable;

public class GenomeLongRegion implements Serializable,Comparable<GenomeLongRegion> {
    private Long start;
    private Long end;

    public GenomeLongRegion(){}
    public GenomeLongRegion(long start,long end){
        this.start=start;
        this.end=end;
    }

    public Long getEnd() {
        return end;
    }

    public void setEnd(Long end) {
        this.end = end;
    }

    public Long getStart() {
        return start;
    }

    public void setStart(Long start) {
        this.start = start;
    }
    public String toString(){
        return start+"\t"+end;
    }

    @Override
    public int compareTo(GenomeLongRegion o) {
        int cop = (int) (this.start - o.getStart());//第一列排序
        if (cop != 0)
            return cop;
        else
            return (int) (this.getEnd()-o.getEnd());
    }
}
