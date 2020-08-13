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
package org.bgi.flexlab.gaea.data.structure.reads;

public class ReadBasicInformation {
	protected String readSequence;
	protected String qualityString;
	protected int MINIMUM_BASE_QUALITY = 33;

	public ReadBasicInformation() {
		this.readSequence = null;
		this.qualityString = null;
	}

	public ReadBasicInformation(String readSequence, String qualityString) {
		this.readSequence = readSequence;
		this.qualityString = qualityString;
	}

	public void setMinimumBaseQuality(int minimumBaseQuality) {
		MINIMUM_BASE_QUALITY = minimumBaseQuality;
	}

	public int getReadLength() {
		return readSequence.length();
	}

	public char getBaseFromRead(int position) {
		return readSequence.charAt(position);
	}

	public String getReadsSequence() {
		return readSequence;
	}

	public char getBaseQuality(int position) {
		return qualityString.charAt(position);
	}

	public int getBaseQualityValue(int position) {
		return (int) getBaseQuality(position) - MINIMUM_BASE_QUALITY;
	}

	public byte getBinaryBase(int i) {
		return (byte) ((readSequence.charAt(i) >> 1) & 0x07);
	}

	public String getQualityString() {
		return qualityString;
	}

	public String getQualityValue() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < qualityString.length(); i++) {
			sb.append(getBaseQualityValue(i));
		}
		return sb.toString();
	}

	public void trim(int leftTrim, int rightTrim) {
		int length = getReadLength();
		if (leftTrim > 0 || rightTrim > 0) {
			if (leftTrim + rightTrim >= length)
				throw new RuntimeException(
						"Trim length is bigger than reads length!");
			
			length -= rightTrim;
			readSequence = readSequence.substring(leftTrim, length);
			qualityString = qualityString.substring(leftTrim, length);
		}
	}
}
