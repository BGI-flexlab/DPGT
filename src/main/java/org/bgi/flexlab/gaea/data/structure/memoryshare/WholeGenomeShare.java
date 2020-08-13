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
package org.bgi.flexlab.gaea.data.structure.memoryshare;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.util.LineReader;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;

public abstract class WholeGenomeShare {
	protected static String DISTRIBUTE_CACHE_FLAG = "distribute.cache.flag";

	public static boolean distributeCache(String chrList, Job job, String cacheName)
			throws IOException, URISyntaxException {
		job.addCacheFile(new URI(chrList + "#" + cacheName));

		Configuration conf = job.getConfiguration();
		Path refPath = new Path(chrList);
		FileSystem fs = refPath.getFileSystem(conf);
		FSDataInputStream refin = fs.open(refPath);
		LineReader in = new LineReader(refin);
		Text line = new Text();

		String chrFile = "";
		String[] chrs = new String[3];
		while ((in.readLine(line)) != 0) {
			chrFile = line.toString();
			chrs = chrFile.split("\t");
			File fileTest = new File(chrs[1]);
			if (fileTest.isFile()) {
				chrs[1] = "file://" + chrs[1];
			}
			job.addCacheFile(new URI(chrs[1] + "#" + chrs[0]));
		}
		in.close();
		refin.close();
		return true;
	}

	protected void loadChromosomeList(String cacheName) {
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(new File(cacheName)));
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e.toString());
		}
		
		try {
			String line;
			while((line = br.readLine()) != null) {
				String[] chrs = line.split("\t");
				// insert chr
				if(!addChromosome(chrs[0])) {
					br.close();
					throw new RuntimeException("map Chromosome "+chrs[1]+" Failed.");
				}
				setChromosome(chrs[0],chrs[0],Integer.parseInt(chrs[2]));
			}
		} catch (NumberFormatException | IOException e) {
			throw new RuntimeException(e.toString());
		}
		try {
			br.close();
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
	
	protected void loadChromosomeList(Path refPath) throws NumberFormatException, IOException{
		Configuration conf = new Configuration();
		FileSystem fs = refPath.getFileSystem(conf);
		FSDataInputStream refin = fs.open(refPath);
		LineReader in = new LineReader(refin);
		Text line = new Text();
		
		String chrFile = "";
		String[] chrs;
		while((in.readLine(line)) != 0){
			chrFile = line.toString();
			chrs = chrFile.split("\t");
			
			// insert chr
			if(!addChromosome(chrs[0])) {
				in.close();
				throw new RuntimeException("map Chromosome "+chrs[1]+" Failed.");
			}
			setChromosome(chrs[1],chrs[0],Integer.parseInt(chrs[2]));
		}
		in.close();
	}

	public static void distributeCacheReference(String chrList, Job job, String cacheName, String distributeCacheFlag) {
		try {
			if (distributeCache(chrList, job, cacheName)) {
				job.getConfiguration().setBoolean(DISTRIBUTE_CACHE_FLAG + "." + cacheName, true);
			}
		} catch (IOException | URISyntaxException e) {
			throw new RuntimeException(e.toString());
		}
	}

	public boolean loadGenome(Configuration conf, String cacheName) {
		boolean isDistributeRef = conf.getBoolean(DISTRIBUTE_CACHE_FLAG + "." + cacheName, false);
		if (!isDistributeRef)
			return false;

		try {
			loadChromosomeList(cacheName);
		} catch (Exception e) {
			throw new RuntimeException(e.toString());
		}
		return true;
	}

	public boolean loadGenome(String refList) {
		try {
			loadChromosomeList(refList);
		} catch (Exception e) {
			return false;
		}
		return true;
	}
	
	public abstract void clean();

	public abstract boolean addChromosome(String chrName);
	
	public abstract void setChromosome(String path,String chrName,int length);
}
