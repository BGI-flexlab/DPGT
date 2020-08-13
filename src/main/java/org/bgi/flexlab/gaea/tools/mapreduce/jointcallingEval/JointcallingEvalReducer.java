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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.VcfLineWritable;

import java.io.IOException;
import java.util.*;

public class JointcallingEvalReducer extends Reducer<Text, VcfLineWritable, NullWritable, Text> {

	private Text resultValue = new Text();
	private List<String> sampleNames;
	private HashMap<String, VCFCodec> vcfCodecs;
	private Path inputPath;
	JointcallingEvalOptions options;

	//	test, baseline, intersection
	private HashMap<String, int[]> stat;

	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		options = new JointcallingEvalOptions();
		Configuration conf = context.getConfiguration();
		options.getOptionsFromHadoopConf(conf);
		sampleNames = new ArrayList<>();
		vcfCodecs = new HashMap<>();
		stat = new HashMap<>();
		stat.put("Total", new int[9]);

		inputPath = new Path(options.getInputFilePath());
		FileSystem fs = inputPath.getFileSystem(conf);
		SingleVCFHeader singleVcfHeader = new SingleVCFHeader();
		singleVcfHeader.readHeaderFrom(inputPath, fs);
		VCFHeader vcfHeader = singleVcfHeader.getHeader();
		VCFHeaderVersion vcfVersion = singleVcfHeader.getVCFVersion(vcfHeader);
		VCFCodec vcfcodec = new VCFCodec();
		vcfcodec.setVCFHeader(vcfHeader, vcfVersion);
		vcfCodecs.put(inputPath.getName(), vcfcodec);
		sampleNames.addAll(vcfHeader.getSampleNamesInOrder());

		Path baselinePath = new Path(options.getBaselineFile());
		fs = baselinePath.getFileSystem(conf);
		singleVcfHeader = new SingleVCFHeader();
		singleVcfHeader.readHeaderFrom(baselinePath, fs);
		vcfHeader = singleVcfHeader.getHeader();
		vcfVersion = singleVcfHeader.getVCFVersion(vcfHeader);
		vcfcodec = new VCFCodec();
		vcfcodec.setVCFHeader(vcfHeader, vcfVersion);
		vcfCodecs.put(baselinePath.getName(), vcfcodec);

		sampleNames.retainAll(vcfHeader.getSampleNamesInOrder());
	}

	@Override
	protected void reduce(Text key, Iterable<VcfLineWritable> values, Context context)
			throws IOException, InterruptedException {
		VariantContext testVariantContext = null;
		VariantContext baselineVariantContext = null;
		Iterator<VcfLineWritable> iter =  values.iterator();
		int count = 0;

		String[] keys = key.toString().split("-");
		String posTag = keys[0] + ":" + keys[1];

		while(iter.hasNext()) {
			VcfLineWritable vcfInput = iter.next();
			String tag = vcfInput.getFileName();
			if (tag.equals(inputPath.getName())) {
				testVariantContext = vcfCodecs.get(tag).decode(vcfInput.getVCFLine());
//				System.out.println(testVariantContext.getReference().toString() + ":" + testVariantContext.getAlleles().toString()+ ":" + testVariantContext.getNAlleles());
				stat.get("Total")[0] += 1;
			} else {
				baselineVariantContext = vcfCodecs.get(tag).decode(vcfInput.getVCFLine());
				stat.get("Total")[1] += 1;
			}
			count++;
		}
		if(count != 2)
			System.out.println("Count is: " + count );

		if(testVariantContext != null && baselineVariantContext != null && testVariantContext.hasSameAllelesAs(baselineVariantContext)){
			stat.get("Total")[2] ++;
		}

		for (String sampleName : sampleNames){
			boolean printDiffPos = false;
			if(stat.containsKey(sampleName)){
				if(testVariantContext != null && baselineVariantContext != null) {
					Genotype testGenotype = testVariantContext.getGenotype(sampleName);
					Genotype baselineGenotype = baselineVariantContext.getGenotype(sampleName);
					if(isVar(testGenotype) && isVar(baselineGenotype)){
						if(sameGenotype(testGenotype, baselineGenotype))
							stat.get(sampleName)[2]++;
						else{
							if(options.isOutputdiff())
								printDiffPos = true;
						}


					}
					if(isVar(testGenotype)) {
						stat.get(sampleName)[0]++;
						if(options.isOutputdiff() && !isVar(baselineGenotype))
							printDiffPos = true;

					}
					if(isVar(baselineGenotype)) {
						stat.get(sampleName)[1]++;
						if(options.isOutputdiff() && !isVar(testGenotype))
							printDiffPos = true;
					}
				}else if(testVariantContext != null){
					Genotype gt = testVariantContext.getGenotype(sampleName);
					if(isVar(gt)) {
						stat.get(sampleName)[0]++;
						if(options.isOutputdiff())
							printDiffPos = true;
					}
				}else if(baselineVariantContext != null){
					Genotype gt = baselineVariantContext.getGenotype(sampleName);
					if(isVar(gt)) {
						stat.get(sampleName)[1]++;
						if(options.isOutputdiff())
							printDiffPos = true;
					}
				}
				if(printDiffPos)
					printDiff(context, sampleName, posTag);
			}else {
				stat.put(sampleName, new int[3]);
			}

		}


	}

	private boolean isVar(Genotype gt){
		if(gt == null || gt.isNoCall() || gt.isHomRef())
			return false;

		for (Allele allele:gt.getAlleles()){
			if(allele.isReference() || Allele.wouldBeStarAllele(allele.getBases()))
				continue;
			return true;
		}
		return false;
	}

	public static boolean sameGenotype(Genotype gt1, Genotype gt2){
		Set<Allele> thisAlleles = new TreeSet<>();
		for (Allele allele: gt1.getAlleles()){
			if(allele.isReference() || Allele.wouldBeStarAllele(allele.getBases()))
				continue;
			thisAlleles.add(allele);
		}
		Set<Allele> otherAlleles = new TreeSet<>();
		for (Allele allele: gt2.getAlleles()){
			if(allele.isReference() || Allele.wouldBeStarAllele(allele.getBases()))
				continue;
			otherAlleles.add(allele);
		}
		return thisAlleles.equals(otherAlleles);
	}

	private void printDiff(Context context, String sampleName, String posTag) throws IOException, InterruptedException {
		resultValue.set("diff\t"+sampleName+"\t"+posTag);
		context.write(NullWritable.get(), resultValue);
	}

	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {
		resultValue.set("REF\t" + stat.get("Total")[0]+"\t"+stat.get("Total")[1] +"\t"+stat.get("Total")[2]);
		context.write(NullWritable.get(), resultValue);
		for (String sampleName : sampleNames) {
			if(!stat.containsKey(sampleName)){
				return;
			}
			resultValue.set(sampleName + "\t" + stat.get(sampleName)[0]+"\t"+stat.get(sampleName)[1] +"\t"+stat.get(sampleName)[2]);
			context.write(NullWritable.get(), resultValue);
		}
	}

}
