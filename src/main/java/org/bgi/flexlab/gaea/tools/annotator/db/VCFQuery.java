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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class VCFQuery extends DBQuery {

	@Override
	public Results query(Condition condition)throws IOException{
		List<String> fields = condition.getFields();
		Results results = new Results();
		List<String> alts = condition.getAlts();
		String key = condition.getConditionString();
		List<HashMap<String,String>> resultList = dbAdapter.getResult( key, fields);

		if (resultList.isEmpty())
			return null;

		for(HashMap<String,String> result :resultList) {
			String resultAltStr = result.get("ALT");
			if (resultAltStr == null) {
				System.err.println("Alt is null:" + condition.getRefTable().getTable() + ". Key:" + key);
				return null;
			}


			String[] resultAlts = resultAltStr.split(",");
			List<HashMap<String, String>> annoResults = splitResult(result, resultAlts.length);
			for (int i = 0; i < resultAlts.length; i++) {
				String alt = resultAlts[i].toUpperCase();
				if (alts.contains(alt)) {
					results.add(alt, annoResults.get(i));
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
