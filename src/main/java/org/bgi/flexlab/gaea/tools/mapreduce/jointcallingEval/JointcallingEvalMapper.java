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
package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingEval;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.InputSplit;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.VcfLineWritable;

import java.io.IOException;
import java.util.HashMap;

public class JointcallingEvalMapper extends Mapper<LongWritable, Text, Text, VcfLineWritable> {

	private Text resultKey;
	private VcfLineWritable resultValue;
	private Configuration conf;
	private HashMap<String, VCFCodec> vcfCodecs;
	private JointcallingEvalOptions options;

	@Override
	protected void setup(Context context)
			throws IOException, InterruptedException {
		resultKey = new Text();
		resultValue = new VcfLineWritable();
		conf = context.getConfiguration();
		vcfCodecs = new HashMap<>();

		options = new JointcallingEvalOptions();
		Configuration conf = context.getConfiguration();
		options.getOptionsFromHadoopConf(conf);

		Path inputPath = new Path(options.getInputFilePath());
		FileSystem fs = inputPath.getFileSystem(conf);
		SingleVCFHeader singleVcfHeader = new SingleVCFHeader();
		singleVcfHeader.readHeaderFrom(inputPath, fs);
		VCFHeader vcfHeader = singleVcfHeader.getHeader();
		VCFHeaderVersion vcfVersion = singleVcfHeader.getVCFVersion(vcfHeader);
		VCFCodec vcfcodec = new VCFCodec();
		vcfcodec.setVCFHeader(vcfHeader, vcfVersion);
		vcfCodecs.put(inputPath.getName(), vcfcodec);

		Path baselinePath = new Path(options.getBaselineFile());
		fs = baselinePath.getFileSystem(conf);
		singleVcfHeader = new SingleVCFHeader();
		singleVcfHeader.readHeaderFrom(baselinePath, fs);
		vcfHeader = singleVcfHeader.getHeader();
		vcfVersion = singleVcfHeader.getVCFVersion(vcfHeader);
		vcfcodec = new VCFCodec();
		vcfcodec.setVCFHeader(vcfHeader, vcfVersion);
		vcfCodecs.put(baselinePath.getName(), vcfcodec);

	}

	private boolean filteVariant(VariantContext variantContext) {
		boolean filter = !variantContext.isVariant() || variantContext.isFiltered();
		if(filter) return true;
		if(variantContext.isSymbolicOrSV())
			return true;
		if(options.getMode().equals(JointcallingEvalOptions.Mode.SNP))
			if(!variantContext.isSNP())
				return true;
		if(options.getMode().equals(JointcallingEvalOptions.Mode.INDEL))
			if(!variantContext.isIndel())
				return true;
		return false;
	}

	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {



		InputSplit inputSplit = context.getInputSplit();
		String fileName = ((FileSplit) inputSplit).getPath().getName();
		System.out.println(fileName);
		VCFCodec vcfcodec = vcfCodecs.get(fileName);
		String vcfLine = value.toString();
		if (vcfLine.startsWith("#")) return;
		VariantContext variantContext = vcfcodec.decode(vcfLine);

		if(filteVariant(variantContext))
			return;

		String chr;

		if(variantContext.getContig().startsWith("chr")){
			chr = variantContext.getContig().substring(3);
		}
		else
			chr = variantContext.getContig();


		resultValue.set(fileName, vcfLine);
//		System.out.println("mapper: "+chr+"-"+variantContext.getStart() + "-" + variantContext.getReference().toString());
		resultKey.set(chr+"-"+variantContext.getStart() + "-" +variantContext.getEnd()+ "-"+ variantContext.getAlternateAllele(0).getBaseString());
//		System.out.println("mapper: " + resultKey.toString() + " " + vcfLine);
		/*根据chr-start-end*/
		context.write(resultKey, resultValue);

	}
	
	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {

	}
}
