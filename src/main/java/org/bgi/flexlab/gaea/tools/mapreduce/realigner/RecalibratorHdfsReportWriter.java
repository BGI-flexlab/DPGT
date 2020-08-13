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
package org.bgi.flexlab.gaea.tools.mapreduce.realigner;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.Path;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFileManager;
import org.bgi.flexlab.gaea.tools.recalibrator.report.RecalibratorReportTable;
import org.bgi.flexlab.gaea.tools.recalibrator.report.RecalibratorReportWriter;

import java.io.IOException;

public class RecalibratorHdfsReportWriter implements RecalibratorReportWriter{
	private FSDataOutputStream stream = null;
	
	public RecalibratorHdfsReportWriter(String out){
		Path path = new Path(out);
		Configuration conf = new Configuration();
		
		stream = HdfsFileManager.getOutputStream(path, conf);
	}

	@Override
	public void write(RecalibratorReportTable table) {
		try {
			stream.write(table.toString().getBytes());
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	@Override
	public void close() {
		try {
			stream.close();
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	@Override
	public void write(String header) {
		try {
			stream.write(header.getBytes());
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
}
