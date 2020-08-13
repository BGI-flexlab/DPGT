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
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.annotator;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.mapreduce.writable.VcfLineWritable;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.io.IOException;
import java.util.HashMap;

public class AnnotationMapper extends Mapper<LongWritable, Text, Text, VcfLineWritable> {

	private Text resultKey;
	private VcfLineWritable resultValue;
	private Configuration conf;
	private HashMap<String, VCFCodec> vcfCodecs;
	@Override
	protected void setup(Context context)
			throws IOException, InterruptedException {
		resultKey = new Text();
		resultValue = new VcfLineWritable();
		conf = context.getConfiguration();
		vcfCodecs = new HashMap<>();

		AnnotatorOptions options = new AnnotatorOptions();
		options.getOptionsFromHadoopConf(conf);

		Path inputPath = new Path(options.getInputFilePath());
		FileSystem fs = inputPath.getFileSystem(conf);
		FileStatus[] files = fs.listStatus(inputPath);
		for(FileStatus file : files) {
			//System.out.println(file.getPath());

			if (file.isFile()) {
				SingleVCFHeader singleVcfHeader = new SingleVCFHeader();
				singleVcfHeader.readHeaderFrom(file.getPath(), fs);
				VCFHeader vcfHeader = singleVcfHeader.getHeader();
				VCFHeaderVersion vcfVersion = singleVcfHeader.getVCFVersion(vcfHeader);
				VCFCodec vcfcodec = new VCFCodec();
				vcfcodec.setVCFHeader(vcfHeader, vcfVersion);
				vcfCodecs.put(file.getPath().getName(), vcfcodec);
			}
		}

//		TODO split vcfLine


	}

	private boolean filteVariant(VariantContext variantContext){
		return !variantContext.isVariant();
	}

	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {
		InputSplit inputSplit = context.getInputSplit();
		String fileName = ((FileSplit) inputSplit).getPath().getName();
		VCFCodec vcfcodec = vcfCodecs.get(fileName);
		String vcfLine = value.toString();
		if (vcfLine.startsWith("#")) return;
		VariantContext variantContext = vcfcodec.decode(vcfLine);

		if(filteVariant(variantContext))
			return;

		String chr = ChromosomeUtils.getNoChrName(variantContext.getContig());

		resultValue.set(fileName, vcfLine);

		int startPrefix = variantContext.getStart()/1000;
		int startRemainder = variantContext.getStart()%1000;
		if(startRemainder <= 5 && startPrefix > 0){
			resultKey.set(chr+"-"+(startPrefix-1));
			context.write(resultKey, resultValue);
		}

		if(startRemainder >= 995 && startPrefix > 0){
			resultKey.set(chr+"-"+(startPrefix+1));
			context.write(resultKey, resultValue);
		}

		resultKey.set(chr+"-"+startPrefix);
//		System.out.println("mapper: " + resultKey.toString() + " " + vcfLine);
		/*根据chr-start-end*/
		context.write(resultKey, resultValue);

	}
	
	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {

	}
}
