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
package org.bgi.flexlab.gaea.tools.annotator.db;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.hadoop.hbase.util.MD5Hash;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo;
import org.bgi.flexlab.gaea.tools.annotator.util.Tuple;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class AnnoDBQuery extends DBQuery {

	public static String INDEX_ALT_COLUMN_NAME = "ALT";

	@Override
	public Results query(Condition condition)throws IOException{
		Results results = new Results();

		HashMap<String,String> result = dbAdapter.getResult(condition.getRefTable().getIndexTable(), condition.getConditionString());
		List<String> alts = condition.getAlts();
		if (result ==null || result.isEmpty()) return null;
		String[] indexAlts = result.get(INDEX_ALT_COLUMN_NAME).split(";");
		for(String alt: alts)
			if(!ArrayUtils.contains(indexAlts, alt))
				return null;

		String mainKeyStr = result.get(condition.getRefTable().getKey());
		String[] mainKeys = mainKeyStr.split(";");

		for(String alt: alts){
			int index = ArrayUtils.indexOf(indexAlts,alt);
			String[] altMainKeys = mainKeys[index].split(",");
			for(String altmk: altMainKeys){
				HashMap<String,String> annoResult = dbAdapter.getResult(condition.getRefTable().getTable(), altmk, condition.getFields());
				if (annoResult ==null || annoResult.isEmpty()){
					System.err.println("Cann't find value from table:"+condition.getRefTable().getTable()+". Key:"+altmk);
					return null;
				}
				results.add(alt, annoResult);
			}
		}
		return results;
	}

	public boolean initInsert(Condition condition,	List<Map<String,String>>
			annos)throws IOException{
		Map<String, String> keys = new HashMap<>();
		int count = 0;
		for(Map<String,String> anno: annos){
			String mainKey = getHashRowKey(condition.getConditionString(), count);
			String alt =  anno.get(INDEX_ALT_COLUMN_NAME);
			dbAdapter.insert(condition.getRefTable().getTable(), mainKey, anno);
			if(keys.containsKey(alt)){
				String mainKeyStr = keys.get(alt) + "," + mainKey;
				keys.put(alt, mainKeyStr);
			}else {
				keys.put(alt, mainKey);
			}
			count ++;
		}

		Map<String, String> indexKV = new HashMap<>();
		Tuple<String, String> tp = Tuple.transMapToTuple(keys);
		indexKV.put(INDEX_ALT_COLUMN_NAME, tp.getFirst());
		indexKV.put(condition.getRefTable().getKey(), tp.getSecond());
		dbAdapter.insert(condition.getRefTable().getIndexTable(), condition.getConditionString(), indexKV);
		return true;
	}

	public boolean addInsert(Condition condition, List<Map<String,String>>
			annos, HashMap<String,String> result)throws IOException{
		if(result == null || result.isEmpty())
			return initInsert(condition, annos);

		Map<String, String> keys = new HashMap<>();
		StringBuilder altStr = new StringBuilder(result.get(INDEX_ALT_COLUMN_NAME));
		String[] alts = altStr.toString().split(";");
		int count = 0;
		for(Map<String,String> anno: annos){

			String alt =  anno.get(INDEX_ALT_COLUMN_NAME);
			int index = ArrayUtils.indexOf(alts, alt);
			if(index != -1){
				continue;
			}

			String mainKey = getHashRowKey(condition.getConditionString(), count + alts.length);
			dbAdapter.insert(condition.getRefTable().getTable(), mainKey, anno);
			if(keys.containsKey(alt)){
				keys.put(alt, keys.get(alt) + "," + mainKey);
			}else {
				keys.put(alt, mainKey);
			}
			count ++;
		}

		StringBuilder mainKeyStr = new StringBuilder(result.get(condition.getRefTable().getKey()));
		Map<String, String> indexKV = new HashMap<>();
		for (String key: keys.keySet()){
			altStr.append(";").append(key);
			mainKeyStr.append(";").append(keys.get(key));
		}
		indexKV.put(INDEX_ALT_COLUMN_NAME, altStr.toString());
		indexKV.put(condition.getRefTable().getKey(), mainKeyStr.toString());
		dbAdapter.insert(condition.getRefTable().getIndexTable(), condition.getConditionString(), indexKV);
		return true;
	}

	public boolean insert(Condition condition,	List<Map<String,String>>
			annos )throws IOException{
		HashMap<String,String> result = dbAdapter.getResult(condition.getRefTable().getIndexTable(), condition.getConditionString());
		return addInsert(condition, annos, result);
	}

	public String getHashRowKey(String key, int count){
		Random random = new Random();
		String mainKey = MD5Hash.getMD5AsHex(Bytes.toBytes(key)).substring(0,10);
		mainKey += String.format("%03x", random.nextInt(2457));
		mainKey += String.format("%03x", count);
		return mainKey;
	}

	@Override
	public void connection(String dbName, DatabaseInfo.DbType dbType, String connInfo) throws IOException{
		dbAdapter = DBAdapterFactory.createDbAdapter(dbType, connInfo);
		dbAdapter.connection("d");
	}

}
