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
package org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.variantrecalibratioin;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.input.bed.RegionHdfsParser;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.*;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import java.io.File;
import java.io.IOException;

public class VariantRecalibrationMapper extends Mapper<LongWritable, VariantContextWritable, IntWritable, Text>{

	private VCFQualityControlOptions options;
	
	private ResourceManager manager;
	
	private GenomeLocationParser genomeLocParser;
	
	private RegionHdfsParser region = null;
	
	/**
	 * 任务初始化设置
	 */
	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		options = new VCFQualityControlOptions();
		options.getOptionsFromHadoopConf(context.getConfiguration());
		
		FastaSequenceFile ref = new FastaSequenceFile(new File(options.getReference()), true);
		genomeLocParser = new GenomeLocationParser(ref.getSequenceDictionary());
		ref.close();

		manager = new ResourceManager(options);
		for(String resource : options.getResources()) {
			System.err.println("resource:"+resource);
			TrainData trainData = new TrainData(options.getReference(), resource);
			if(trainData.isDB()) {
				trainData.setType(new DBResource());
			} else {
				trainData.setType(new FileResource());
			}
			trainData.initialize();
			manager.addTrainingSet(trainData);
		}
		if( !manager.checkHasTrainingSet() ) {
			throw new UserException.CommandLineException( "No training set found! Please provide sets of known polymorphic loci marked with the training=true ROD binding tag. For example, -resource:hg19-hapmap" );
		}
		if( !manager.checkHasTruthSet() ) {
			throw new UserException.CommandLineException( "No truth set found! Please provide sets of known polymorphic loci marked with the truth=true ROD binding tag. For example, -resource:hg19-hapmap" );
		}
		
		if(options.getRegion() != null){
			region = new RegionHdfsParser();
            region.parseBedFileFromHDFS(options.getRegion(), false);
		}
	}
	
	@Override
	public void map(LongWritable key, VariantContextWritable value, Context context) throws IOException, InterruptedException {
		VariantContext vc = value.get();
		if(!validContext(vc))
			return;
		
		boolean inRegion = false;
		for(int i = vc.getStart() ; i <= vc.getEnd() ; i++){
			if(region != null && region.isPositionInRegion(vc.getContig(), i - 1)) {
                inRegion = true;
                break;
            }
		}
		
		if(region != null && !inRegion)
			return;
		
		VariantDatumMessenger datum = new VariantDatumMessenger.Builder(manager, vc, options)
														 .decodeAnnotations()
														 .setLoc(genomeLocParser)
														 .setOriginalQual()
														 .setFlagV()
														 .setPrior()
														 .build();
		if(datum != null) {
			context.write(new IntWritable((int)key.get()), new Text(datum.toString()));
		}
	}

	public boolean validContext(VariantContext vc) {
		return vc != null && 
				(vc.isNotFiltered() || options.getIgnoreInputFilters().containsAll(vc.getFilters())) &&
				ResourceManager.checkVariationClass(vc, options.getMode());
	}
	
}
