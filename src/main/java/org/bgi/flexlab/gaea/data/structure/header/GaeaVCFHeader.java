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
package org.bgi.flexlab.gaea.data.structure.header;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

public class GaeaVCFHeader implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4160938012339601002L;

	public static final String VCF_HEADER_PROPERTY = "vcfHeader";
	
	public void writeHeaderToHDFS(String outputPath, Configuration conf){
		writeToHDFS(new Path(conf.get(VCF_HEADER_PROPERTY)));
	}
	
	public boolean writeToHDFS(Path path){
		ObjectOutputStream ostream = null;
		try {
			FileSystem fs = path.getFileSystem(new Configuration());
			ostream = new ObjectOutputStream(fs.create(path));
			ostream.writeObject(this);
			ostream.close();
			fs.close();
		} catch (IOException e) {
			ostream = null;
			e.printStackTrace();
			return false;
		}
		return true;
	}

	public static GaeaVCFHeader readFromFile(String file) throws IOException, ClassNotFoundException {
		Configuration conf = new Configuration();
		Path path = new Path(file);
		FileSystem fs = path.getFileSystem(conf);
		FSDataInputStream stream = fs.open(path);
		ObjectInputStream ostream = new ObjectInputStream(stream);
		GaeaVCFHeader vcfHeader = (GaeaVCFHeader) ostream.readObject();
		ostream.close();
		return vcfHeader;
	}

	public static GaeaVCFHeader loadVcfHeader(boolean cache, Configuration conf){
		GaeaVCFHeader vcfHeaderTmp = null;
		try {
			if(cache){//distribute cache
				vcfHeaderTmp = readFromFile("VcfHeaderObj");
			} else {//no distribute cache
				vcfHeaderTmp = readFromFile(conf.get(VCF_HEADER_PROPERTY));
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return vcfHeaderTmp;
	}
}
