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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.util;

import org.bgi.flexlab.gaea.data.exception.UserException;

import java.util.ArrayList;
import java.util.List;

public class StringUtils {
	public static String join(String[] str, String seperator) {
		boolean first = true;

		StringBuilder sb = new StringBuilder();
		for (String s : str) {
			if (first) {
				sb.append(s);
				first = false;
			} else {
				sb.append(seperator + s);
			}
		}
		return sb.toString();
	}

	public static List<Integer> getWordStarts(String line) {
		if (line == null)
			throw new UserException("line is null");
		List<Integer> starts = new ArrayList<Integer>();
		int stop = line.length();
		for (int i = 1; i < stop; i++)
			if (Character.isWhitespace(line.charAt(i - 1)))
				if (!Character.isWhitespace(line.charAt(i)))
					starts.add(i);
		return starts;
	}

	public static String[] splitFixedWidth(String line) {
		List<Integer> starts = getWordStarts(line);
		
		return splitFixedWidth(line,starts);
	}
	
	public static String[] splitFixedWidth(String line,List<Integer> starts) {
		if (line == null)
			throw new UserException("line is null");
		if (starts == null)
			throw new UserException("columnStarts is null");
		int startCount = starts.size();
		String[] row = new String[startCount + 1];
		if (startCount == 0) {
			row[0] = line.trim();
		} else {
			row[0] = line.substring(0, starts.get(0)).trim();
			for (int i = 1; i < startCount; i++)
				row[i] = line.substring(starts.get(i - 1), starts.get(i)).trim();
			row[startCount] = line.substring(starts.get(startCount - 1)).trim();
		}
		return row;
	}
}
