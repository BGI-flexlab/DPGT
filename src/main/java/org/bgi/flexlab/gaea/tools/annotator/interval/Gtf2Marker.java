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
 * Copyright (C)  2016  Pablo Cingolani(pcingola@users.sourceforge.net)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator.interval;

import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * An interval intended as a mark
 *
 * @author pcingola
 *
 */
public class Gtf2Marker extends GffMarker {

	private static final long serialVersionUID = 5416962964309837838L;

	static final String ATTRIBUTE_PATTERN_REGEX = "\\s*(\\S+)\\s+\"(.*?)\"\\s*;";
	static final Pattern ATTRIBUTE_PATTERN = Pattern.compile(ATTRIBUTE_PATTERN_REGEX);

	public Gtf2Marker() {
		super();
	}

	public Gtf2Marker(Genome genome, String line) {
		super(genome, line);
	}

	public Gtf2Marker(Marker parent, int start, int end, boolean strandMinus, String id) {
		super(parent, start, end, strandMinus, id);
	}

	/**
	 * Parse attributes
	 */
	@Override
	protected void parseAttributes(String attrStr) {
		keyValues = new HashMap<>();
		keys = new HashSet<String>();

		if (attrStr.length() > 0) {
			Matcher matcher = ATTRIBUTE_PATTERN.matcher(attrStr);
			while (matcher.find()) {
				if (matcher.groupCount() >= 2) {
					String key = matcher.group(1).toLowerCase();
					String value = matcher.group(2);
					add(key, value);
				}
			}
		}
	}

}
