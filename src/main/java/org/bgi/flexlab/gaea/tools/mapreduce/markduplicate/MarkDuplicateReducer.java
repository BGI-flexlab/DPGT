package org.bgi.flexlab.gaea.tools.mapreduce.markduplicate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.DuplicationKeyWritable;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.markduplicate.MarkDuplicatesFunc;

import java.io.IOException;
import java.util.ArrayList;

public class MarkDuplicateReducer extends Reducer<DuplicationKeyWritable, SamRecordWritable, NullWritable, SamRecordWritable>{
    private MarkDuplicatesFunc mark = new MarkDuplicatesFunc();
    private SAMFileHeader samHeader;

    @Override
    public void setup(Context context){
        Configuration conf = context.getConfiguration();
        samHeader = SamHdfsFileHeader.getHeader(conf);
    }

    public void reduce(DuplicationKeyWritable key, Iterable<SamRecordWritable> values, Context context) throws IOException, InterruptedException {
        //deal unmapped reads
        if(key.getChrIndex() == -1) {
            for(SamRecordWritable s : values) {
            	GaeaSamRecord sam = new GaeaSamRecord(samHeader,s.get());
                SamRecordWritable w = new SamRecordWritable();
                w.set(sam);
                context.write(NullWritable.get(), w);
            }
            return;
        }

        //collect reads cluster and mark duplicate
        ArrayList<SAMRecord> sams = new ArrayList<>();
        int n = 0;
        for(SamRecordWritable s : values) {
        	GaeaSamRecord sam = new GaeaSamRecord(samHeader,s.get());
            if (n>100000) {
                sam.setDuplicateReadFlag(true);
                SamRecordWritable w = new SamRecordWritable();
                w.set(sam);
                context.write(NullWritable.get(), w);
            }else {
                sams.add(sam);
            }
            n++;
        }

        if(sams.size() > 1)
            mark.markDup(sams);
        for(SAMRecord sam : sams) {
            SamRecordWritable w = new SamRecordWritable();
            w.set(sam);
            context.write(NullWritable.get(), w);
        }
    }
}
