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
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.data.mapreduce.writable.VcfLineWritable;
import org.bgi.flexlab.gaea.data.structure.header.SingleVCFHeader;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.annotator.AnnotationEngine;
import org.bgi.flexlab.gaea.tools.annotator.SampleAnnotationContext;
import org.bgi.flexlab.gaea.tools.annotator.VcfAnnoContext;
import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.db.DBAnnotator;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.io.IOException;
import java.util.*;

public class AnnotationReducer extends Reducer<Text, VcfLineWritable, Text, Text> {

	private Text resultKey = new Text();
	private Text resultValue = new Text();
	private AnnotatorOptions options;
	private HashMap<String, VCFCodec> vcfCodecs;
	private HashMap<String, VCFHeader> vcfHeaders;
	private AnnotationEngine annoEngine;
	private DBAnnotator dbAnnotator;
	ReferenceShare genomeShare;
	Config userConfig;
	long mapTime = 0;
	long mapCount = 0;
	
	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		options = new AnnotatorOptions();
		Configuration conf = context.getConfiguration();
		options.getOptionsFromHadoopConf(conf);

		long start = System.currentTimeMillis();
		genomeShare = new ReferenceShare();
		genomeShare.loadChromosomeList(options.getReferenceSequencePath());
		if(options.isDebug())
			System.err.println("genomeShare耗时：" + (System.currentTimeMillis()-start)+"毫秒");

		userConfig = new Config(conf, genomeShare);

		start = System.currentTimeMillis();
		AnnotatorBuild annoBuild = new AnnotatorBuild(userConfig);
		userConfig.setSnpEffectPredictor(annoBuild.createSnpEffPredictor());
		annoBuild.buildForest();
		if(options.isDebug())
			System.err.println("build SnpEffectPredictor耗时：" + (System.currentTimeMillis()-start)+"毫秒");
		vcfCodecs = new HashMap<>();
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
				vcfCodecs.put(file.getPath().getName(), vcfcodec);
				vcfHeaders.put(file.getPath().getName(), vcfHeader);
			}

		}
		if(options.isDebug())
			System.err.println("getVCFHeader耗时：" + (System.currentTimeMillis()-start)+"毫秒");

		annoEngine = new AnnotationEngine(userConfig);

		start = System.currentTimeMillis();
		//用于从数据库中查找信息
		dbAnnotator = new DBAnnotator(userConfig);
		try {
			dbAnnotator.connection();
		} catch (InstantiationException | IllegalAccessException
				| ClassNotFoundException e) {
			e.printStackTrace();
		}
		System.err.println("dbAnnotator.connection耗时：" + (System.currentTimeMillis()-start)+"毫秒");

	}

	@Override
	protected void reduce(Text key, Iterable<VcfLineWritable> values, Context context)
			throws IOException, InterruptedException {
		long start = System.currentTimeMillis();
		Iterator<VcfLineWritable> iter =  values.iterator();
		Map<Long, VcfAnnoContext> posVariantInfo = new HashMap<>();
		TreeSet<Long> positions = new TreeSet<>();

		while(iter.hasNext()) {
			VcfLineWritable vcfInput =  iter.next();
			String fileName = vcfInput.getFileName();
			String vcfLine = vcfInput.getVCFLine();
			VariantContext variantContext =  vcfCodecs.get(fileName).decode(vcfLine);
			long posKey = (long)variantContext.getStart();
			posKey = (posKey << 32) | variantContext.getEnd();

			if(positions.contains(posKey)){
				posVariantInfo.get(posKey).add(variantContext, fileName);
			}else {
				VcfAnnoContext vcfAnnoContext = new VcfAnnoContext(variantContext, fileName);
				posVariantInfo.put(posKey, vcfAnnoContext);
				positions.add(posKey);
			}
		}

		for (long posKey : positions) {
			VcfAnnoContext vcfAnnoContext = posVariantInfo.get(posKey);
			String chr = ChromosomeUtils.getNoChrName(vcfAnnoContext.getContig());
			String posPrefix = chr + "-" + vcfAnnoContext.getStart() / 1000;

			long higherPositionKey = positions.higher(posKey) != null ? positions.higher(posKey) : (long) 0;
			// 标记附近有其他变异的点
			if (higherPositionKey != 0 && (higherPositionKey >> 32) - vcfAnnoContext.getEnd() <= 5) {
				VcfAnnoContext vcfAnnoContextNear = posVariantInfo.get(higherPositionKey);
				for (SampleAnnotationContext sac : vcfAnnoContext.getSampleAnnoContexts().values()) {
					String sampleName = sac.getSampleName();
					if (vcfAnnoContextNear.hasSample(sampleName)) {
						sac.setHasNearVar();
						vcfAnnoContextNear.getSampleAnnoContexts().get(sampleName).setHasNearVar();
					}
				}
			}

			if (!posPrefix.equals(key.toString()))
				continue;

			if (userConfig.getFields().contains("FLKSEQ")) {
				ChromosomeInformationShare chrShare = genomeShare.getChromosomeInfo(vcfAnnoContext.getContig());
				int lelfStart = vcfAnnoContext.getStart() - 11;
				lelfStart = lelfStart < 0 ? 0 : lelfStart;
				int lelfEnd = vcfAnnoContext.getStart() - 2;
				lelfEnd = lelfEnd < 0 ? 0 : lelfEnd;
				String prefixSeq = chrShare.getGA4GHBaseSequence(lelfStart, lelfEnd);
				int rightStart = vcfAnnoContext.getStart() + vcfAnnoContext.getRefStr().length() - 1;
				String suffixSeq = chrShare.getGA4GHBaseSequence(rightStart, rightStart + 9);
				String flankSeq = prefixSeq + "." + suffixSeq;
				vcfAnnoContext.setFlankSeq(flankSeq);
			}

			if (!options.isUseDatabaseCache() || !dbAnnotator.annotate(vcfAnnoContext, "ANNO")) {
				if (!annoEngine.annotate(vcfAnnoContext)) {
					continue;
				}
				dbAnnotator.annotate(vcfAnnoContext);
				if (options.isDatabaseCache())
					dbAnnotator.insert(vcfAnnoContext, "ANNO");
			}

			if (options.getOutputFormat() == AnnotatorOptions.OutputFormat.VCF) {
				Map<String, List<VariantContext>> annos = vcfAnnoContext.toAnnotationVariantContexts(userConfig.getFieldsWithoutVariant());
				for (String filename : annos.keySet()) {
					resultKey.set(filename);
					VCFHeader vcfHeader = vcfHeaders.get(filename);
					List<VariantContext> variantContexts = annos.get(filename);
					VCFEncoder vcfEncoder = new VCFEncoder(vcfHeader, true, true);
					for (VariantContext vc : variantContexts) {
						String vcfLine = vcfEncoder.encode(vc);
						resultValue.set(vcfLine);
						context.write(resultKey, resultValue);
					}
				}
			} else {
				List<String> annoLines = vcfAnnoContext.toAnnotationStrings(userConfig.getFields());
				for (String annoLine : annoLines) {
					resultValue.set(annoLine);
					context.write(resultKey, resultValue);
				}
			}
		}

		if(options.isDebug()) {
			System.err.println("step3:" + (System.currentTimeMillis() - start) + "ms");
			mapTime += System.currentTimeMillis() - start;
			mapCount++;
		}
	}

	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {
		dbAnnotator.disconnection();
		if(options.isDebug())
			System.err.println("dbAnnotator平均耗时(mapTime/mapCount)：" +mapTime+"/"+mapCount+" = ? 毫秒");
	}
}
