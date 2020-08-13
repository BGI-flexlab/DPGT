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
package org.bgi.flexlab.gaea.data.structure.dbsnp;

import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.Job;
import org.bgi.flexlab.gaea.data.structure.memoryshare.WholeGenomeShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class DbsnpShare extends WholeGenomeShare {
	private static final String CACHE_NAME = "dbsnpList";

	private Map<String, ChromosomeDbsnpShare> dbsnpInfo = new ConcurrentHashMap<String, ChromosomeDbsnpShare>();
	
	public DbsnpShare(String dbsnpPath,String refPath){
		indexExist(dbsnpPath,refPath);
	}
	
	public static void indexExist(String dbsnpPath,String refPath){
		String chrList = dbsnpPath + VcfIndex.INDEX_SUFFIX;
		if(chrList.startsWith("file://"))
			chrList = chrList.substring("file://".length());
		File file = new File(chrList);
		
		if(!file.exists()){
			VcfIndex index = new VcfIndex();
			index.buildIndex(refPath, dbsnpPath, null);
		}
	}

	public static boolean distributeCache(String chrList, Job job) {
		try {
			return distributeCache(chrList, job, CACHE_NAME);
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		} catch (URISyntaxException e) {
			throw new RuntimeException(e.toString());
		}
	}

	public void loadChromosomeList() {
		loadChromosomeList(CACHE_NAME);
	}

	public void loadChromosomeList(String chrList) {
		try {
			chrList = "file://"+chrList;
			loadChromosomeList(new Path(chrList));
		} catch (IllegalArgumentException | IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	@Override
	public boolean addChromosome(String chrName) {
		ChromosomeDbsnpShare newChr = new ChromosomeDbsnpShare();
		dbsnpInfo.put(chrName, newChr);
		if (dbsnpInfo.get(chrName) != null) {
			return true;
		} else {
			return false;
		}
	}

	public ChromosomeDbsnpShare getChromosomeDbsnp(String chrName) {
		chrName = ChromosomeUtils.formatChrName(chrName);
		if (!dbsnpInfo.containsKey(chrName))
			return null;
		return dbsnpInfo.get(chrName);
	}

	/*public long getStartPosition(String chrName, int winNum, int winSize) {
		chrName = ChromosomeUtils.formatChrName(chrName);
		if (!dbsnpInfo.containsKey(chrName)){
			return -1;
		}

		return dbsnpInfo.get(chrName).getStartPosition(winNum, winSize);
	}*/
	
	public long getStartPosition(String chrName, int winNum, int winSize) {
		return getStartPosition(chrName,winNum,winNum+1,winSize);
	}
	
	public long getStartPosition(String chrName, int startWinNum,int endWinNum, int winSize) {
		chrName = ChromosomeUtils.formatChrName(chrName);
		if (!dbsnpInfo.containsKey(chrName)){
			return -1;
		}
		return dbsnpInfo.get(chrName).getStartPosition(startWinNum, endWinNum,winSize);
	}

	public long getStartPosition(String chrName, int winNum) {
		return getStartPosition(chrName, winNum, 0);
	}

	public Map<String, ChromosomeDbsnpShare> getDbsnpMap() {
		return this.dbsnpInfo;
	}

	public void setChromosome(String path, String chrName, int length) {
		if (dbsnpInfo.containsKey(chrName)) {
			// map chr and get length
			dbsnpInfo.get(chrName).loadChromosome(path);
			dbsnpInfo.get(chrName).setLength(length);
			dbsnpInfo.get(chrName).setChromosomeName(chrName);
		}
	}

	@Override
	public void clean() {
		for(String key : dbsnpInfo.keySet()){
			ChromosomeDbsnpShare share = dbsnpInfo.get(key);
			try {
				share.clean();
			} catch (Exception e) {
				throw new RuntimeException(e.toString());
			}
		}
	}
}
