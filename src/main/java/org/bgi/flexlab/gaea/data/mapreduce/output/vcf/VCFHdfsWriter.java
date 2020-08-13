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
package org.bgi.flexlab.gaea.data.mapreduce.output.vcf;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFFileWriter;

import java.io.BufferedOutputStream;
import java.io.IOException;

public class VCFHdfsWriter extends VCFFileWriter{
	/**
	 * 
	 */
	private static final long serialVersionUID = 4557833627590736290L;

	public VCFHdfsWriter(String filePath, boolean doNotWriteGenotypes,
			final boolean allowMissingFieldsInHeader, Configuration conf) throws IOException {
		super(filePath, doNotWriteGenotypes, allowMissingFieldsInHeader, conf);
	}
	
	@Override
	public void initOutputStream(String filePath, Configuration conf) {
		os = new BufferedOutputStream(HdfsFileManager.getOutputStream(new Path(filePath), conf));
	}
}
