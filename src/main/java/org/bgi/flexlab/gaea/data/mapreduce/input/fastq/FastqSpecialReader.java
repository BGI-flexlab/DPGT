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
package org.bgi.flexlab.gaea.data.mapreduce.input.fastq;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.FileSplit;

import java.io.IOException;

public class FastqSpecialReader extends FastqBasicReader {
	private String readFlag = "1";

	public FastqSpecialReader(Configuration job, FileSplit split,
			byte[] recordDelimiter) throws IOException {
		super(job, split, recordDelimiter);
		readFlag = getReadFlagFromFileName(split.getPath().getName());
	}

	private String getReadFlagFromFileName(String fileName) {
		String flag = "1";
		int index = fileName.lastIndexOf(".fq");
		String tempFlag = fileName.substring(index - 1, index);
		if (tempFlag.equals("1") || tempFlag.equals("2"))
			flag = tempFlag;
		return flag;
	}

	@Override
	public boolean next(Text key, Text value) throws IOException {
		if (key == null) {
			key = new Text();
		}
		if (value == null) {
			value = new Text();
		}
		int newSize = 0;
		boolean iswrongFq = false;
		while (pos < end) {
			Text tmp = new Text();
			String[] st = new String[4];
			int startIndex = 0;

			if (firstLine != "") {
				st[0] = firstLine;
				st[1] = secondLine;
				startIndex = 2;
				firstLine = "";
				secondLine = "";
			}

			for (int i = startIndex; i < 4; i++) {
				newSize = in.readLine(tmp, maxLineLength, Math.max(
						(int) Math.min(Integer.MAX_VALUE, end - pos),
						maxLineLength));

				if (newSize == 0) {
					iswrongFq = true;
					break;
				}
				pos += newSize;
				st[i] = tmp.toString();
			}
			if (!iswrongFq) {
				st[0] = st[0] + "/" + readFlag;
				int index = st[0].lastIndexOf("/");
				String tempkey = st[0].substring(1, index).trim();
				if(sampleID != null && sampleID != "+")
					st[2] = sampleID;
				
				key.set(tempkey);
				value.set(st[0] + "\t" + st[1] + "\t" + st[2] + "\t" + st[3]);
			} else {
				LOG.warn("wrong fastq reads:blank line among fq file or end of file!");
			}
			break;
		}
		if (newSize == 0 || iswrongFq) {
			key = null;
			value = null;
			return false;
		} else {
			return true;
		}
	}

}
