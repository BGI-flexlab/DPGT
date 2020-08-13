package org.bgi.flexlab.gaea.data.mapreduce.output.vcf;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.RecordWriter;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.seqdoop.hadoop_bam.LazyBCFGenotypesContext;
import org.seqdoop.hadoop_bam.LazyParsingGenotypesContext;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext;
import org.seqdoop.hadoop_bam.VCFRecordWriter;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public abstract class GaeaVCFRecordWriter<K> extends RecordWriter<K,VariantContextWritable> {

	private VCFCodec codec = new VCFCodec();
	private VariantContextWriter writer;
	private boolean writerHeaderIsInit = false;
	private VCFHeader header;

	private LazyVCFGenotypesContext.HeaderDataCache vcfHeaderDataCache =
		new LazyVCFGenotypesContext.HeaderDataCache();
	private LazyBCFGenotypesContext.HeaderDataCache bcfHeaderDataCache =
		new LazyBCFGenotypesContext.HeaderDataCache();

	/** A VCFHeader is read from the input Path. */
	public GaeaVCFRecordWriter(
			Path output, Path input, boolean writeHeader, TaskAttemptContext ctx)
		throws IOException
	{
		final AsciiLineReader r = new AsciiLineReader(
			input.getFileSystem(ctx.getConfiguration()).open(input));

		final Object h = codec.readHeader(new AsciiLineReaderIterator(r));
		if (!(h instanceof VCFHeader))
			throw new IOException("No VCF header found in "+ input);

		r.close();

		init(output, (VCFHeader)h, writeHeader, ctx);
	}
	public GaeaVCFRecordWriter(
			Path output, VCFHeader header, boolean writeHeader,
			TaskAttemptContext ctx)
		throws IOException
	{
		init(
			output.getFileSystem(ctx.getConfiguration()).create(output),
			header, writeHeader);
	}
	public GaeaVCFRecordWriter(
			OutputStream output, VCFHeader header, boolean writeHeader)
		throws IOException
	{
		init(output, header, writeHeader);
	}

	// Working around not being able to call a constructor other than as the
	// first statement...
	private void init(
			Path output, VCFHeader header, boolean writeHeader,
			TaskAttemptContext ctx)
		throws IOException
	{
		init(
			output.getFileSystem(ctx.getConfiguration()).create(output),
			header, writeHeader);
	}
	private void init(
			OutputStream output, VCFHeader header, boolean writeHeader)
		throws IOException
	{
		final StoppableOutputStream stopOut =
			new StoppableOutputStream(!writeHeader, output);

		writer = VariantContextWriterFactory.create(
			stopOut, null, VariantContextWriterFactory.NO_OPTIONS);

		stopOut.stopped = false;

		if(header != null)
			setInputHeader(this.header = header);
	}
	
	private void initWriter(VCFHeader vcfHeader){
		if(!writerHeaderIsInit){
			writer.writeHeader(vcfHeader);
			writerHeaderIsInit = true;
		}
	}

	@Override 
	public void close(TaskAttemptContext ctx) throws IOException {
		writer.close();
	}

	/** Used for lazy decoding of genotype data. Of course, each input record
	 * may have a different header, but we currently only support one header
	 * here... This is in part due to the fact that it's not clear what the best
	 * solution is. */
	public void setInputHeader(VCFHeader header) {
		vcfHeaderDataCache.setHeader(header);
		bcfHeaderDataCache.setHeader(header);
	}

	protected void writeRecord(VariantContext vc) {
		final GenotypesContext gc = vc.getGenotypes();
		if (gc instanceof LazyParsingGenotypesContext)
			((LazyParsingGenotypesContext)gc).getParser().setHeaderDataCache(
				gc instanceof LazyVCFGenotypesContext ? vcfHeaderDataCache
				                                      : bcfHeaderDataCache);

		writer.add(vc);
	}
	
	protected void writeRecord(VariantContext vc,VCFHeader header) {
		initWriter(header);
		setInputHeader(header);
		writeRecord(vc);
	}
}

final class StoppableOutputStream extends FilterOutputStream {
	public boolean stopped;

	public StoppableOutputStream(boolean startStopped, OutputStream out) {
		super(out);
		stopped = startStopped;
	}

	@Override public void write(int b) throws IOException {
		if (!stopped) super.write(b);
	}
	@Override public void write(byte[] b) throws IOException {
		if (!stopped) super.write(b);
	}
	@Override public void write(byte[] b, int off, int len) throws IOException {
		if (!stopped) super.write(b, off, len);
	}
}
