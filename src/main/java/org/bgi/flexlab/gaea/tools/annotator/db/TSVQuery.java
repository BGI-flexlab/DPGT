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

import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;

public class TSVQuery extends DBQuery {

	private static final long serialVersionUID = 805441802476012341L;

	/**
	 * 根据condition查询数据库
	 * @param condition
	 * @return
	 * @throws IOException
	 */
	public Results query(Condition condition)throws IOException{
		List<String> fields = condition.getFields();
		Results results = new Results();

		HashMap<String,String> result = dbAdapter.getResult(condition.getRefTable().getIndexTable(), condition.getConditionString());
		if (result ==null || result.isEmpty()) return null;
		List<String> alts = condition.getAlts();

		String keyStr = result.get(condition.getConditionString());

		String[] keys = keyStr.split(",");
		for (String key : keys) {
			result = dbAdapter.getResult("data", key);

			HashMap<String,String> annoResult = new HashMap<>();
			for (String field : fields) {
				annoResult.put(field, result.get(field));
			}

			if (result ==null || result.isEmpty()){
				System.err.println("Cann't find value from table:"+condition.getRefTable().getTable()+". Key:"+key);
				return null;
			}

			String resultAltStr = result.get("ALT");
			if (resultAltStr == null) {
				System.err.println("Alt is null:"+condition.getRefTable().getTable()+". Key:"+key);
				return null;
			}

			if (!resultAltStr.contains(",")) {
				resultAltStr = resultAltStr.toUpperCase();
				if(alts.contains(resultAltStr)){
					results.add(resultAltStr, annoResult);
				}
			}else {
				String[] resultAlts = resultAltStr.split(",");
				List<HashMap<String, String>> annoResults = splitResult(annoResult, resultAlts.length);
				for (int i = 0; i < resultAlts.length; i++) {
					String alt = resultAlts[i].toUpperCase();
					if(alts.contains(alt)){
						results.add(alt, annoResults.get(i));
					}
				}
			}
		}

		return results;
	}

	public void connection(String dbName, DatabaseInfo.DbType dbType, String connInfo) throws IOException{
		dbAdapter = DBAdapterFactory.createDbAdapter(dbType, connInfo);
		dbAdapter.connection(dbName);
	}
	
}
