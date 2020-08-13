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

import org.bgi.flexlab.gaea.tools.annotator.VcfAnnoContext;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo.DbType;
import org.bgi.flexlab.gaea.tools.annotator.AnnotationContext;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;

public class DBQuery implements Serializable {

	private static final long serialVersionUID = -897843908487603204L;
	
	DBAdapter dbAdapter = null;
	Results results = null;
	Condition condition = null;
	
	/**
	 * 执行query并判断结果，通过getResults方法获取结果
	 * @param condition
	 * @return
	 * @throws IOException
	 */
	boolean executeQuery(Condition condition) throws IOException {
		this.condition = condition;
		Results results = query(condition);
		
		if (results == null || results.isEmpty()) {
			return false;
		}else {
			this.results = results;
			return true;
		}
	}
	
	/**
	 * 根据condition查询数据库
	 * @param condition
	 * @return
	 * @throws IOException
	 */
	public Results query(Condition condition)throws IOException{
		Results results = new Results();

		HashMap<String,String> result = dbAdapter.getResult(condition.getRefTable().getIndexTable(), condition.getConditionString());
		if (result ==null || result.isEmpty()) return null;
		List<String> alts = condition.getAlts();
		String keyStr = result.get(condition.getRefTable().getKey());
		String[] keys = keyStr.split(",");
		for (String key : keys) {
			result = dbAdapter.getResult(condition.getRefTable().getTable(), key);

			if (result ==null || result.isEmpty())
				return null;

			String resultAltStr = result.get("ALT");
			if (resultAltStr == null) {
				System.err.println("Alt is null:"+condition.getRefTable().getTable()+". Key:"+key);
				return null;
			}

			HashMap<String,String> annoResult = new HashMap<>();
			for(String field: condition.getFields()){
				annoResult.put(field, result.get(field));
			}

			String[] resultAlts = resultAltStr.split(",");
			List<HashMap<String, String>> annoResults = splitResult(annoResult, resultAlts.length);
			for (int i = 0; i < resultAlts.length; i++) {
				String alt = resultAlts[i].toUpperCase();
				if(alts.contains(alt)){
					results.add(alt, annoResults.get(i));
				}
			}
		}
		
		return results;
	}
	
	/**
	 * 对含多个变异的结果进行分割
	 */
	protected List<HashMap<String, String>> splitResult(HashMap<String, String> result, int altNum) {
		List<HashMap<String, String>> resultList = new ArrayList<>();
		if(altNum == 1) {
			resultList.add(result);
			return resultList;
		}

		for (int i = 0; i < altNum; i++) {
			resultList.add(new HashMap<>());
		}

		for (Entry<String, String> entry : result.entrySet()) {
			String v = entry.getValue();
			if(v == null) continue;
			if(v.contains(",")){
				String[] values = v.split(",");
				if(altNum == values.length){
					for (int i = 0; i < altNum; i++) {
						resultList.get(i).put(entry.getKey(), values[i]);
					}
				}else {
					for (int i = 0; i < altNum; i++) {
						resultList.get(i).put(entry.getKey(), v);
					}
				}
			}else {
				for (int i = 0; i < altNum; i++) {
					resultList.get(i).put(entry.getKey(), entry.getValue());
				}
			}
		}
		return resultList;
	}
	
	public LinkedList<HashMap<String,String>> getResultList(String tag) {
		return results.get(tag);
	}

	public HashMap<String,String> getMergeResult(AnnotationContext ac) {
		return results.getMergeResult(ac.getAlt());
	}

	public boolean insert(Condition condition,	List<Map<String,String>> annos)throws IOException{
		return true;
	}

	/**
	 * 对查询结果results进行矫正
	 */
	//abstract void adjustResult(HashMap<String,String> result);

	Results getResults(){
		return results;
	}

	public void disconnection() throws IOException {
		dbAdapter.disconnection();
	}

	public void connection(String dbName, DbType dbType, String connInfo) throws IOException{
		dbAdapter = DBAdapterFactory.createDbAdapter(dbType, connInfo);
		dbAdapter.connection("data");
	}
}
