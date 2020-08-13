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

import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformation;
import org.bgi.flexlab.gaea.data.structure.vcf.AbstractVCFLoader.PositionalVariantContext;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.io.*;
import java.util.Arrays;

public class VcfIndex extends ReferenceIndex {
	public final static String INDEX_SUFFIX = ".window.idx";
	public final static int WINDOW_SIZE = 100;
	private int capacity = Long.SIZE / Byte.SIZE;

	public VcfIndex() {
	}

	private void saveAsBinary(String outputPath, byte[] binaries) throws IOException {
		System.err.println("starting writing:"+binaries.length);
		DataOutputStream out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputPath)));
		out.write(binaries);
		out.close();
	}

	private void setBinaryIndex(int winNum, byte[] indexs, long position) {
		int index = winNum * capacity;

		int move;
		for (int i = 0; i < capacity; i++) {
			move = 56 - i * 8;
			indexs[index + i] = (byte) ((position >> move) & 0xff);
		}
	}

	private void insertSnpInformation(ChromosomeInformation curChrInfo, VariantContext context, byte[] sequences) {
		for (int pos = context.getStart(); pos <= context.getEnd(); pos++) {
			if (pos > curChrInfo.getLength())
				break;
			curChrInfo.insertSnpInformation(pos - 1);
		}
	}

	private void fileWriter(FileWriter bnListWriter, String outputPath, String lastChrName, int lastLength) {
		try {
			bnListWriter.write(lastChrName);
			bnListWriter.write("\t");
			bnListWriter.write(outputPath + "/" + lastChrName + ".dbsnp.bn");
			bnListWriter.write("\t");
			bnListWriter.write(String.valueOf(lastLength));
			bnListWriter.write("\n");
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	@Override
	protected void dbsnpParser(String dbsnpPath, String outputPath) {
		if(dbsnpPath == null)
			return;
		VCFLocalLoader reader = null;
		try {
			reader = new VCFLocalLoader(dbsnpPath);
		} catch (IOException e2) {
			throw new RuntimeException(e2.toString());
		}

		String lastChrName = null;
		int lastLength = -1;
		ChromosomeInformation curChrInfo = null;
		byte[] binaryIndex = null;

		String dbsnpList = dbsnpPath + INDEX_SUFFIX;
		FileWriter bnListWriter = null;
		try {
			bnListWriter = new FileWriter(new File(dbsnpList));
		} catch (IOException e1) {
			throw new RuntimeException(e1.toString());
		}

		int lastWinNum = -1;

		while (reader.hasNext()) {
			PositionalVariantContext posContext = reader.next();
			long currPos = posContext.getPosition();
			VariantContext context = posContext.getVariantContext();
			String chrName = ChromosomeUtils.formatChrName(context.getContig());

			if (lastChrName == null || !lastChrName.equals(chrName)) {
				curChrInfo = chromosomeInfoMap.get(chrName);
				if (null == curChrInfo)
					throw new RuntimeException(
							"> Failed appending dbSNP information. No information related to chromosome name: "
									+ chrName);
				int len = curChrInfo.getLength() / WINDOW_SIZE;
				if ((curChrInfo.getLength() % WINDOW_SIZE) != 0)
					len++;

				if (lastChrName != null && !lastChrName.equals(chrName)) {
					try {
						saveAsBinary(outputPath + "/" + lastChrName + ".dbsnp.bn", binaryIndex);
					} catch (IOException e) {
						throw new RuntimeException(e.toString());
					}
					fileWriter(bnListWriter, outputPath, lastChrName, lastLength);
				}
				
				binaryIndex = null;
				binaryIndex = new byte[len * capacity];
				Arrays.fill(binaryIndex, 0, binaryIndex.length, (byte) 0);

				lastChrName = chrName;
				lastLength = curChrInfo.getLength();
			}

			int currWinNum = ((context.getStart() - 1) / WINDOW_SIZE);

			if (lastWinNum != currWinNum) {
				lastWinNum = currWinNum;
				setBinaryIndex(currWinNum, binaryIndex, currPos);
			}

			insertSnpInformation(curChrInfo, context, binaryIndex);
		}

		try {
			saveAsBinary(outputPath + "/" + lastChrName + ".dbsnp.bn", binaryIndex);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}

		fileWriter(bnListWriter, outputPath, lastChrName, lastLength);

		reader.close();

		try {
			bnListWriter.close();
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	public static void main(String[] args) {
		VcfIndex index = new VcfIndex();

		if (args.length < 2) {
			System.err.println("java -Xmx10g -jar gaea-1.0.0.jar reference_path dbsnp_path output_path");
			System.exit(1);
		}
		
		if(args.length == 2)
			index.buildIndex(args[0], args[1], null);
		else
			index.buildIndex(args[0], args[1], args[2]);

		System.out.println("build index finish!!");
	}
}
