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

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.util.SAMUtils;

import java.util.Arrays;

public class ReadBasicCompressionInformation {
	protected byte[] readBases;
	protected byte[] qualities;
	//protected byte[] compressedReadBases;
	protected static int MINIMUM_BASE_QUALITY = 33;

	public ReadBasicCompressionInformation() {
		this.readBases = null;
		this.qualities = null;
	}

	public ReadBasicCompressionInformation(ReadBasicCompressionInformation read) {
		readBases = Arrays.copyOf(read.readBases, read.readBases.length);
		qualities = Arrays.copyOf(read.qualities, read.qualities.length);
	}

	public ReadBasicCompressionInformation(byte[] readBases, byte[] qualities) {
		this.readBases = readBases;
		this.qualities = qualities;
	}
	
	public ReadBasicCompressionInformation(String readSequence, String qualityString) {
		this.readBases = SAMUtils.bytesToCompressedBasesGaea(readSequence.getBytes());
		this.qualities = qualityString.getBytes();
	}

	/**
	 *  set minimum base quality
	 * @param minimumBaseQuality
	 */
	public void setMinimumBaseQuality(int minimumBaseQuality) {
		MINIMUM_BASE_QUALITY = minimumBaseQuality;
	}

	/**
	 * get read length
	 * @return
	 */
	public int getReadLength() {
		return qualities.length;
	}

	/**
	 * get base
	 * @param position
	 * @return
	 */
	public byte getReadBase(int position) {
		return readBases[position];
	}

	/**
	 * get base index
	 * @param position start from 0
	 * @return
	 */
	public byte getBinaryBase(int position) {
		return (byte) ((getReadBase(position) & 0x0f) >> 1);
	}

	/**
	 * get compressed read bases
	 * @return
	 */
	public byte[] getReadBases() {
		return readBases;
	}


	/**
	 *  get read bases in String
	 * @return read bases in String
	 */
	public String getReadBasesString() {
		return String.valueOf(getReadBases());
	}

	/**
	 * set read bases
	 * @param readBases
	 */
	public void setReadBases(byte[] readBases) {
		this.readBases = readBases;
	}

	/**
	 * get read region seq
	 * @param qStart
	 * @param length
	 * @return
	 */
	public char[] getReadBases(int qStart, int length) {
		char[] bases = new char[length];

		for(int i = 0; i < length; i++) {
			bases[i] = (char)readBases[qStart + i];
		}
		return bases;
	}

	public String getReadBasesString(int qstart, int length) {
		return String.valueOf(getReadBases(qstart, length));
	}

	/**
	 * compress read bases.
	 * @return
	 */
	public byte[] getCompressedReadBasesBytes() {
		//return SAMUtils.compressedBasesToBytesGaea(qualities.length, readBases, 0);
		return SAMUtils.bytesToCompressedBasesGaea(readBases);
	}

	public void getReadBasesFromCompressedReadBasesBytes(byte[] compressedReadBases, int readLength) {
		readBases = SAMUtils.compressedBasesToBytesGaea(readLength, compressedReadBases, 0);
	}

	/**
	 * get qualities bytes
	 * @return
	 */
	public byte[] getQualities() {
		return qualities;
	}

	/**
	 * set qualities bytes
	 * @param qualities
	 */
	public void setQualities(byte[] qualities) {
		this.qualities = qualities;
	}

	public String getQualitiesString() {
		return String.valueOf(qualities);
	}

	/**
	 * get base quality at position
	 * @param position
	 * @return
	 */
	public byte getBaseQuality(int position) {
		return qualities[position];
	}

	/**
	 * get base quality true value
	 * @param position
	 * @return
	 */
	public byte getBaseQualityValue(int position) {
		return (byte) (getBaseQuality(position) - MINIMUM_BASE_QUALITY);
	}

	public boolean isEmpty() {
		return readBases == null || qualities.length == 0;
	}

	public void emptyRead() {
		readBases = null;
		qualities = null;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		return sb.toString();
	}
}
