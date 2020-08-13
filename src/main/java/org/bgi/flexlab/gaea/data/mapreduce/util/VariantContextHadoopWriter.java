package org.bgi.flexlab.gaea.data.mapreduce.util;

import java.io.IOException;

import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.output.vcf.GaeaVariantContextWriter;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextHadoopWriter implements GaeaVariantContextWriter{

	@SuppressWarnings("rawtypes")
	private Context context = null;
	
	private VariantContextWritable writable = new VariantContextWritable();
	
	private VCFHeader header = new VCFHeader();
	
	public VariantContextHadoopWriter() {}
	
	@SuppressWarnings("rawtypes")
	public VariantContextHadoopWriter(Context context,VCFHeader header) {
		this.context = context;
		this.header = header;
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public void write(VariantContext record) {
		if(record != null) {
			writable.set(record, header);
			try {
				context.write(NullWritable.get(), writable);
			} catch (IOException e) {
				throw new UserException(e.toString());
			} catch (InterruptedException e) {
				throw new UserException(e.toString());
			}
		}
	}
	
	@Override
	public void close() {
	}

	@Override
	public void add(VariantContext vc) {
		write(vc);
	}

	@SuppressWarnings("rawtypes")
	public Context getContext() {
		return this.context;
	}
}
