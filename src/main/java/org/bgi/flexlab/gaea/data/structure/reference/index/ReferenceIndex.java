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
package org.bgi.flexlab.gaea.data.structure.reference.index;

import org.bgi.flexlab.gaea.data.exception.NullFilePathException;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformation;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.io.*;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

/**
 * create reference and dbSNP file index
 *
 * @author ZhangYong, ZhangZhi
 */
public abstract class ReferenceIndex {
	private final static String HDFS_PRIX = "file://"; 
	/**
	 * chromosome information map
	 */
	protected Map<String, ChromosomeInformation> chromosomeInfoMap = new ConcurrentHashMap<String, ChromosomeInformation>();

	/**
	 * read fasta reference
	 * 
	 * @param refSeqPath
	 *            reference file path
	 * @throws IOException
	 */
	protected void loadReference(String refSeqPath) throws IOException {
		if (null != refSeqPath && !refSeqPath.equals("")) {
			File in = new File(refSeqPath);
			BufferedReader reader = new BufferedReader(new FileReader(in));

			String line = null;
			StringBuilder refSeq = new StringBuilder();
			String curChrName = "";
			// 读取参考基因组序列
			while ((line = reader.readLine()) != null) {
				// '>'字符是一个染色体序列信息的起始标志
				if (line.length() != 0 && '>' == line.charAt(0)) {
					// 处理之前的一个染色体
					if (chromosomeInfoMap.containsKey(curChrName)) {
						// 二进制化染色体碱基序列
						chromosomeInfoMap.get(curChrName).setBinarySequence(refSeq.toString());
						System.out.println("> Finished loading chromosome: " + curChrName);
					}

					// 获取染色体名称
					int pos;
					for (pos = 1; pos != line.length() && '\t' != line.charAt(pos) && '\n' != line.charAt(pos)
							&& ' ' != line.charAt(pos) && '\r' != line.charAt(pos) && '\f' != line.charAt(pos); pos++) {
					}
					String newChrName = ChromosomeUtils.formatChrName(line.substring(1, pos));

					// 判断添加染色体信息是否成功
					if (!addChromosome(newChrName)) {
						StringBuilder errorDescription = new StringBuilder();
						errorDescription.append("> Insert Chromosome ");
						errorDescription.append(newChrName);
						errorDescription.append(" Failed.");
						System.err.println(errorDescription.toString());
					}
					curChrName = newChrName;
					refSeq.setLength(0);
				} else {
					refSeq.append(line);
				}
			}
			// 处理最后一个染色体的信息
			if (refSeq.length() != 0 && chromosomeInfoMap.containsKey(curChrName)) {
				// 二进制化染色体碱基序列
				chromosomeInfoMap.get(curChrName).setBinarySequence(refSeq.toString());
			}
			reader.close();
		} else {
			throw new NullFilePathException("input", "reference");
		}
	}

	/**
	 * get chromosome information for chrName
	 */
	protected ChromosomeInformation getChromosomeInformation(String chrName) {
		return chromosomeInfoMap.get(chrName);
	}

	/**
	 * 向chromosomeInfoMap添加染色体信息
	 */
	protected boolean addChromosome(String chrName) {
		ChromosomeInformation newChr = new ChromosomeInformation();
		chromosomeInfoMap.put(chrName, newChr);
		return (chromosomeInfoMap.get(chrName) != null) ? true : false;
	}

	protected void referenceSaveAsBinary(String refIndexOutputPath) throws IOException {
		StringBuilder outputRefListPath = new StringBuilder(); // 编码后的二进制库文件的索引文件路径
		outputRefListPath.append(refIndexOutputPath);
		outputRefListPath.append("/ref_bn.list");
		FileWriter refList = new FileWriter(new File(outputRefListPath.toString()));
		Iterator<Entry<String, ChromosomeInformation>> iter = chromosomeInfoMap.entrySet().iterator();
		while (iter.hasNext()) {
			Entry<String, ChromosomeInformation> entry = iter.next();
			String chrName = entry.getKey();
			ChromosomeInformation curChrInfo = entry.getValue();
			StringBuilder outputFileName = new StringBuilder();
			outputFileName.append(refIndexOutputPath);
			outputFileName.append("/");
			outputFileName.append(chrName);
			outputFileName.append(".fa.bn");
			curChrInfo.outputChrInformation(outputFileName.toString()); // 输出编码后的二进制库文件

			int length = curChrInfo.getLength();
			refList.write(chrName);
			refList.write("\t");
			refList.write(outputFileName.toString());
			refList.write("\t");
			refList.write(String.valueOf(length));
			refList.write("\n");
		}
		refList.close();
	}
	
	private String absoluteOutputPath(String output){
		File file = new File(output);
		String abs = file.getAbsolutePath();
		return abs;
	}
	
	private void createDirectory(String path){
		File file = new File(path);
		
		if(!file.exists())
			file.mkdirs();
		else
			throw new RuntimeException("output directory is exist,please remove first!");
	}

	public void buildIndex(String refPath, String dbsnpListPath, String indexOutputPath) {
		if(refPath.startsWith(HDFS_PRIX))
			refPath = refPath.substring(HDFS_PRIX.length());
		
		if(dbsnpListPath.startsWith(HDFS_PRIX))
			dbsnpListPath = dbsnpListPath.substring(HDFS_PRIX.length());
		
		if(indexOutputPath == null){
			File file = new File(refPath);
			indexOutputPath = file.getParentFile()+"/WindowIndex/";
		}
			
		indexOutputPath = absoluteOutputPath(indexOutputPath);
		if (!indexOutputPath.endsWith("/"))
			indexOutputPath += "/";
		
		String referenceIndexOutputPath = indexOutputPath + "reference";
		String dbsnpIndexOutputPath = indexOutputPath + "dbsnp";
		
		createDirectory(referenceIndexOutputPath);
		createDirectory(dbsnpIndexOutputPath);

		try {
			loadReference(refPath);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}

		if(dbsnpListPath != null)
			dbsnpParser(dbsnpListPath,dbsnpIndexOutputPath);

		try {
			referenceSaveAsBinary(referenceIndexOutputPath);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
	
	protected abstract void dbsnpParser(String dbSnpListPath,String outputPath);
}
