package org.bgi.flexlab.gaea.data.mapreduce.output.vcf;

import java.io.IOException;
import java.io.OutputStream;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.seqdoop.hadoop_bam.VariantContextWithHeader;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GaeaKeyIgnoringVCFRecordWriter<K> extends GaeaVCFRecordWriter<K> {
	
	public GaeaKeyIgnoringVCFRecordWriter(
			Path output, Path input, boolean writeHeader, TaskAttemptContext ctx)
		throws IOException{
		super(output, input, writeHeader, ctx);
	}
	
	public GaeaKeyIgnoringVCFRecordWriter(
			Path output, VCFHeader header, boolean writeHeader,
			TaskAttemptContext ctx)
		throws IOException{
		super(output, header, writeHeader, ctx);
	}
	
	public GaeaKeyIgnoringVCFRecordWriter(
			OutputStream output, VCFHeader header, boolean writeHeader)
		throws IOException{
		super(output, header, writeHeader);
	}

	@Override 
	public void write(K ignored, VariantContextWritable vc) {
		VariantContext v = vc.get();
		if(v instanceof VariantContextWithHeader)
			writeRecord(v,((VariantContextWithHeader)v).getHeader());
		else
			writeRecord(v);
	}
}
