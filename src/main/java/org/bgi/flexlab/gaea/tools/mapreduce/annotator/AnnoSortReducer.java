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

import htsjdk.variant.vcf.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.bgi.flexlab.gaea.data.mapreduce.writable.PairWritable;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.tools.annotator.config.Config;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class AnnoSortReducer extends Reducer<PairWritable, Text, NullWritable, Text> {

	private MultipleOutputs<NullWritable,Text> multipleOutputs = null;
	private Text resultValue;
	private boolean printHeader = true;
	private AnnotatorOptions options;
	private HashMap<String, VCFHeader> vcfHeaders;
	private String annoFieldNameHeader;

	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		resultValue = new Text();
		multipleOutputs = new MultipleOutputs<>(context);
		Configuration conf = context.getConfiguration();
		Config userConfig = new Config(conf);
		options = new AnnotatorOptions();
		options.getOptionsFromHadoopConf(conf);
		annoFieldNameHeader = userConfig.getVCFHeaderString();

		vcfHeaders = new HashMap<>();
		Path inputPath = new Path(options.getInputFilePath());
		FileSystem fs = inputPath.getFileSystem(conf);
		FileStatus[] files = fs.listStatus(inputPath);

		for(FileStatus file : files) {
			if (file.isFile()) {
				SingleVCFHeader singleVcfHeader = new SingleVCFHeader();
				singleVcfHeader.readHeaderFrom(file.getPath(), fs);
				VCFHeader vcfHeader = singleVcfHeader.getHeader();
				VCFHeaderVersion vcfVersion = SingleVCFHeader.getVCFHeaderVersion(vcfHeader);
				VCFCodec vcfcodec = new VCFCodec();
				vcfcodec.setVCFHeader(vcfHeader, vcfVersion);
				vcfHeaders.put(file.getPath().getName(), vcfHeader);
			}
		}
	}

	@Override
	protected void reduce(PairWritable key, Iterable<Text> values, Context context)
			throws IOException, InterruptedException {
		//TODO fix like AnnotationSortReducer
		if(printHeader) {
			VCFHeader vcfHeader = vcfHeaders.get(key.getFirst());
			VCFHeaderLine annoVcfHeaderLine = new VCFInfoHeaderLine("ANNO", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "fieldName:"+annoFieldNameHeader);
			vcfHeader.addMetaDataLine(annoVcfHeaderLine);
			for(VCFHeaderLine vcfHeaderLine: vcfHeader.getMetaDataInInputOrder()){
				resultValue.set(VCFHeader.METADATA_INDICATOR + vcfHeaderLine.toString());
				multipleOutputs.write(NullWritable.get(), resultValue, key.getFirst());
			}

			StringBuilder fieldLine = new StringBuilder(VCFHeader.HEADER_INDICATOR);
			boolean isFirst = true;
			for (final VCFHeader.HEADER_FIELDS field : vcfHeader.getHeaderFields() ) {
				if ( isFirst )
					isFirst = false; // don't write out a field separator
				else
					fieldLine.append(VCFConstants.FIELD_SEPARATOR);
				fieldLine.append(field.toString());
			}

			if ( vcfHeader.hasGenotypingData() ) {
				fieldLine.append(VCFConstants.FIELD_SEPARATOR);
				fieldLine.append("FORMAT");
				for (final String sample : vcfHeader.getGenotypeSamples() ) {
					fieldLine.append(VCFConstants.FIELD_SEPARATOR);
					fieldLine.append(sample);
				}
			}

			resultValue.set(fieldLine.toString());
			multipleOutputs.write(NullWritable.get(), resultValue, key.getFirst());
			printHeader = false;
		}
		for (Text inputLine : values) {
			resultValue.set(inputLine);
			multipleOutputs.write(NullWritable.get(), resultValue, key.getFirst());
		}
	}

	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {
		multipleOutputs.close();
	}
}
