package org.bgi.flexlab.gaea.tools.mapreduce.jointcalling;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.RecordReader;
import org.apache.hadoop.mapreduce.TaskAttemptContext;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.seqdoop.hadoop_bam.LazyParsingGenotypesContext;
import org.seqdoop.hadoop_bam.LazyVCFGenotypesContext;
import org.seqdoop.hadoop_bam.VariantContextWritable;
import org.seqdoop.hadoop_bam.util.MurmurHash3;

import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class JointCallingVCFRecordReader extends RecordReader<LongWritable, VariantContextWritable> {
	private final LongWritable key = new LongWritable();
	private final VariantContextWritable vc = new VariantContextWritable();

	private VCFCodec codec = new VCFCodec();
	private AsciiLineReaderIterator it;
	private AsciiLineReader reader;

	private long length;

	private final Map<String, Integer> contigDict = new HashMap<String, Integer>();
	
	private VCFHeader header = null;
	
	private LazyVCFGenotypesContext.HeaderDataCache vcfHeaderDataCache =
    		new LazyVCFGenotypesContext.HeaderDataCache();

	@Override
	public void initialize(InputSplit spl, TaskAttemptContext ctx) throws IOException {
		final FileSplit split = (FileSplit) spl;

		this.length = split.getLength();

		final Path file = split.getPath();
		final FileSystem fs = file.getFileSystem(ctx.getConfiguration());

		final FSDataInputStream ins = fs.open(file);

		reader = new AsciiLineReader(ins);
		it = new AsciiLineReaderIterator(reader);

		final Object h = codec.readHeader(it);
		if (!(h instanceof FeatureCodecHeader) || !(((FeatureCodecHeader) h).getHeaderValue() instanceof VCFHeader))
			throw new IOException("No VCF header found in " + file);

		header = (VCFHeader) ((FeatureCodecHeader) h).getHeaderValue();
		
		vcfHeaderDataCache.setHeader(header);

		contigDict.clear();
		for (final VCFContigHeaderLine contig : header.getContigLines())
			contigDict.put(contig.getID(), contig.getContigIndex());

		// Note that we create a new reader here, so reader.getPosition() is 0
		// at
		// start regardless of the value of start. Hence getProgress() and
		// nextKeyValue() don't need to use start at all.
		final long start = split.getStart();
		if (start != 0) {
			ins.seek(start - 1);
			reader = new AsciiLineReader(ins);
			reader.readLine(); // NOTE: skip incomplete line!
			it = new AsciiLineReaderIterator(reader);
		} else { // it seems that newer versions of the reader peek ahead one
					// more line from the input
			long current_pos = it.getPosition();
			ins.seek(0);
			reader = new AsciiLineReader(ins);
			it = new AsciiLineReaderIterator(reader);
			while (it.hasNext() && it.getPosition() <= current_pos && it.peek().startsWith("#")) {
				it.next();
			}
			if (!it.hasNext() || it.getPosition() > current_pos)
				throw new IOException("Empty VCF file " + file);
		}
	}

	@Override
	public void close() {
		reader.close();
	}

	@Override
	public float getProgress() {
		return length == 0 ? 1 : (float) reader.getPosition() / length;
	}

	@Override
	public LongWritable getCurrentKey() {
		return key;
	}

	@Override
	public VariantContextWritable getCurrentValue() {
		return vc;
	}

	@Override
	public boolean nextKeyValue() throws IOException {
		if (!it.hasNext() || it.getPosition() >= length)
			return false;

		final String line = it.next();
		final VariantContext v = codec.decode(line);
		
		VariantContext vtemp = codec.decode(line);
		
		GenotypesContext gc = vtemp.getGenotypes();
		if (gc instanceof LazyParsingGenotypesContext){
			((LazyParsingGenotypesContext)gc).getParser().setHeaderDataCache(vcfHeaderDataCache);
		}
		
		CommonInfo info = v.getCommonInfo();
		if(!info.hasAttribute("SM"))
			info.putAttribute("SM",vtemp.getSampleNamesOrderedByName().get(0));

		Integer chromIdx = contigDict.get(v.getContig());
		if (chromIdx == null)
			chromIdx = (int) MurmurHash3.murmurhash3(v.getContig(), 0);

		key.set((long) chromIdx << 32 | (long) (v.getStart() - 1));
		vc.set(v);
		return true;
	}
}
