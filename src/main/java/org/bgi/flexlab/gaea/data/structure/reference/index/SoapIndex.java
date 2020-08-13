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

import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformation;

import java.io.*;

public class SoapIndex extends ReferenceIndex {

	private static byte getSnpBasicInfo(String[] info) {
		byte snpBasicInfo = 0;

		// IsValidated信息
		if (info[3].trim().equals("1")) {
			snpBasicInfo |= 1;
		}

		// IsHapMap信息
		if (info[2].trim().equals("1")) {
			snpBasicInfo |= 2;
		}

		// FreqG
		if (!info[8].trim().equals("0")) {
			snpBasicInfo |= 4;
		}
		// FreqT
		if (!info[7].trim().equals("0")) {
			snpBasicInfo |= 8;
		}
		// FreqC
		if (!info[6].trim().equals("0")) {
			snpBasicInfo |= 16;
		}
		// FreqA
		if (!info[5].trim().equals("0")) {
			snpBasicInfo |= 32;
		}

		return snpBasicInfo;
	}

	@SuppressWarnings("resource")
	private int soapParser(String dbSNPPath, String outPath, String outIndexPath) {
		int dbsnpNum = 0;
		if (null != dbSNPPath && !dbSNPPath.equals("")) {
			DataOutputStream out = null;
			DataOutputStream index = null;
			File in = new File(dbSNPPath);
			BufferedReader reader = null;
			try {
				out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outPath)));
				index = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outIndexPath)));
				reader = new BufferedReader(new FileReader(in));
			} catch (FileNotFoundException e) {
				throw new RuntimeException(e.toString());
			}

			int pos = -1;
			String currentName = null;
			String[] splitArray = null;
			String line = null;
			// 按行读取dbSNP数据文件
			try {
				while ((line = reader.readLine()) != null && line.length() != 0) {
					// dbSNP数据格式:
					// Chr\tPos\thapmap?\tvalidated?\tis_indel?\tA\tC\tT\tG\trsID\n
					splitArray = line.split("\t");

					// 如果该位点为Indel，则不添加该SNP信息
					if (splitArray[4].trim().equals("1")) {
						continue;
					}

					// 等位基因频率数组
					float[] freq = { Float.parseFloat(splitArray[5].trim()), Float.parseFloat(splitArray[6].trim()),
							Float.parseFloat(splitArray[7].trim()), Float.parseFloat(splitArray[8].trim()) };
					byte alleleCount = 0;
					byte zeroCount = 0;
					int nonZeroPos = 0; // 非零等位基因位置
					// 计算存在于dbSNP中的等位基因个数
					for (int base = 0; base < 4; base++) {
						if (freq[base] > 0) { // 等位基因频率大于0
							alleleCount += 1;
							if (alleleCount == 1)
								nonZeroPos = base;
						} else {
							zeroCount += 1;
						}
					}
					// 存在于dbSNP中的等位基因个数小于2，或等位基因频率都为0时，则不添加该SNP信息
					if (alleleCount < 2 || zeroCount == 4) {
						continue;
					}

					// 染色体名称
					currentName = splitArray[0].trim().toLowerCase();
					if (!currentName.contains("chr")) {
						currentName = "chr" + currentName;
					}

					// 位点坐标
					pos = Integer.parseInt(splitArray[1].trim());
					index.writeInt(pos - 1);
					dbsnpNum++;

					// 设置SNP基本信息
					out.writeByte(getSnpBasicInfo(splitArray));

					// 设置等位基因频率
					out.writeFloat(freq[nonZeroPos]);

					ChromosomeInformation curChrInfo = chromosomeInfoMap.get(currentName);
					if (null != curChrInfo) {
						// 向染色体信息中添加dbSNP信息
						if (pos < curChrInfo.getLength()) {
							pos -= 1;
							curChrInfo.insertSnpInformation(pos);
						}
					} else {
						throw new RuntimeException(
								"> Failed appending dbSNP information. No information related to chromosome name: "
										+ currentName);
					}
				}
			} catch (NumberFormatException e) {
				throw new RuntimeException(e.toString());
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
			try {
				reader.close();
				out.close();
				index.close();
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}

		} else {
			System.err.println("> Empty path of dbSNP file.");
		}
		return dbsnpNum;
	}

	@SuppressWarnings("resource")
	@Override
	protected void dbsnpParser(String dbSnpListPath,String outputPath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(dbSnpListPath));
			String line;

			FileWriter dbsnpIndexWriter = new FileWriter(new File(outputPath + "/dbsnp_bn.list"));
			try {
				while ((line = br.readLine()) != null) {
					String[] lineArray = line.split("\t");
					String outPath = outputPath + "/" + lineArray[0] + ".dbsnp.bn";
					String outIndexPath = outputPath + "/" + lineArray[0] + ".dbsnp.index";
					int dbsnpNumber = soapParser(lineArray[1], outPath, outIndexPath);

					dbsnpIndexWriter.write(lineArray[0]);
					dbsnpIndexWriter.write("\t");
					dbsnpIndexWriter.write(outPath);
					dbsnpIndexWriter.write("\t");
					dbsnpIndexWriter.write(outIndexPath);
					dbsnpIndexWriter.write("\t");
					dbsnpIndexWriter.write(String.valueOf(dbsnpNumber));
					dbsnpIndexWriter.write("\n");
				}
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
			
			br.close();
			dbsnpIndexWriter.close();
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e.toString());
		} catch (IOException e1) {
			throw new RuntimeException(e1.toString());
		}
	}
}
