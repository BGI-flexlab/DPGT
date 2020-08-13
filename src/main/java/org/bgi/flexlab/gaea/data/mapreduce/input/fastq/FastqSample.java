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

import org.apache.hadoop.fs.Path;

public class FastqSample {
	private int id;
	private String fileName;
	private String sampleName;
	private String readsGroup;
	private String fastq1;
	private String fastq2;
	private String adapter1;
	private String adapter2;
	private String indexString;
	
	private void setDefault() {
		fastq1 = null;
		fastq2 = null;
		adapter1 = null;
		adapter2 = null;
		indexString = null;
	}

	public FastqSample() {
		setDefault();
	}

	public boolean setSampleList(String line, boolean keepFile) {
		String splitArray[] = line.split("\t");

		if (splitArray.length < 6) {
			return false;
		}

		id = Integer.parseInt(splitArray[0]);
		readsGroup = splitArray[1];
		sampleName = getRGtag(readsGroup, "SM");
		fileName = getReportFileName(readsGroup);

		// format path
		splitArray[2] = new Path(splitArray[2]).toString();
		splitArray[3] = new Path(splitArray[3]).toString();
		if (!keepFile) {
			for (int i = 2; i < 6; i++) {
				if (splitArray[i].startsWith("file:///")) {
					splitArray[i] = splitArray[i].substring(7);
				} else {
					if (splitArray[i].startsWith("file:/")) {
						splitArray[i] = splitArray[i].substring(5);
					}
				}
			}
		}
		if (!splitArray[2].toLowerCase().equals("null"))
			fastq1 = splitArray[2];
		if (!splitArray[3].toLowerCase().equals("null"))
			fastq2 = splitArray[3];
		if (!splitArray[4].toLowerCase().equals("null"))
			adapter1 = splitArray[4];
		if (!splitArray[5].toLowerCase().equals("null"))
			adapter2 = splitArray[5];

		if (splitArray.length > 6)
			indexString = splitArray[6];
		return true;
	}

	private String getReportFileName(String rg) {
		StringBuilder fileNamebuilBuilder = new StringBuilder();
		fileNamebuilBuilder.append(sampleName);
		fileNamebuilBuilder.append(".");
		fileNamebuilBuilder.append(getRGtag(rg, "ID"));
		return fileNamebuilBuilder.toString();
	}

	private String getRGtag(String rg, String tag) {
		String temp;
		int index = rg.indexOf(tag + ":");
		if (index >= 0) {
			temp = rg.substring(index);
			index = temp.indexOf("\\t");
			if (index < 0) {
				temp = temp.trim();
				index = temp.length();
			}
			temp = temp.substring(0, index);
			String[] splite = temp.split(":");
			return splite[1];
		} else
			return null;
	}

	/**
	 * @return the readsGroup
	 */
	public String getReadsGroup() {
		return readsGroup;
	}

	/**
	 * @param readsGroup
	 *            the readsGroup to set
	 */
	public void setReadsGroup(String readsGroup) {
		this.readsGroup = readsGroup;
	}

	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	/**
	 * @param id
	 *            the id to set
	 */
	public void setId(int id) {
		this.id = id;
	}

	/**
	 * @return the sampleName
	 */
	public String getSampleName() {
		return sampleName;
	}

	/**
	 * @param sampleName
	 *            the sampleName to set
	 */
	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}

	/**
	 * @return the fastq1
	 */
	public String getFastq1() {
		return fastq1;
	}

	/**
	 * @param fastq1
	 *            the fastq1 to set
	 */
	public void setFastq1(String fastq1) {
		this.fastq1 = fastq1;
	}

	/**
	 * @return the fastq2
	 */
	public String getFastq2() {
		return fastq2;
	}

	/**
	 * @param fastq2
	 *            the fastq2 to set
	 */
	public void setFastq2(String fastq2) {
		this.fastq2 = fastq2;
	}

	/**
	 * @return the adapter1
	 */
	public String getAdapter1() {
		return adapter1;
	}

	/**
	 * @param adapter1
	 *            the adapter1 to set
	 */
	public void setAdapter1(String adapter1) {
		this.adapter1 = adapter1;
	}

	/**
	 * @return the adapter2
	 */
	public String getAdapter2() {
		return adapter2;
	}

	/**
	 * @param adapter2
	 *            the adapter2 to set
	 */
	public void setAdapter2(String adapter2) {
		this.adapter2 = adapter2;
	}

	/**
	 * @return the fileName
	 */
	public String getFileName() {
		return fileName;
	}

	/**
	 * @param fileName
	 *            the fileName to set
	 */
	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public String getIndex() {
		return indexString;
	}
}
