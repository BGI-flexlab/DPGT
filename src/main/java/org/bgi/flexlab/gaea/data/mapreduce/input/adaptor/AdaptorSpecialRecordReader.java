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
package org.bgi.flexlab.gaea.data.mapreduce.input.adaptor;

import org.apache.hadoop.io.Text;

import java.io.IOException;

public class AdaptorSpecialRecordReader extends AdaptorRecordReader {
	@Override
	public boolean nextKeyValue() throws IOException {
		Text temp = new Text();
		while (pos < end) {
			int newSize = in.readLine(temp, maxLineLength,
					Math.max((int) Math.min(Integer.MAX_VALUE, end - pos),
							maxLineLength));
			if (newSize == 0) {
				return false;
			}
			pos += newSize;
			if (newSize < maxLineLength) {
				String line = temp.toString();
				String[] splitLines = line.split("\t");
				if (splitLines[0].startsWith("#")) {
					continue;
				}
				String tempkey, tempvalue;
				splitLines[0] += "/1";
				int index = splitLines[0].lastIndexOf("/");

				tempkey = splitLines[0].substring(0, index).trim();
				tempvalue = splitLines[0].substring(index + 1).trim();

				key.set(tempkey);
				value.set(tempvalue);
				return true;
			}

			// line too long. try again
			LOG.info("Skipped line of size " + newSize + " at pos "
					+ (pos - newSize));
		}
		return false;
	}
}
