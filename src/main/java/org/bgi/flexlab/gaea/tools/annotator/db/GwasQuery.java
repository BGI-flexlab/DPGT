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

import java.io.IOException;
import java.util.HashMap;

public class GwasQuery extends DBQuery {

	private static final long serialVersionUID = 805441802476032672L;

	@Override
	public Results query(Condition condition)throws IOException{
		Results results = new Results();

		HashMap<String,String> result = dbAdapter.getResult(condition.getRefTable().getIndexTable(), condition.getConditionString());
		if (result ==null || result.isEmpty()) return null;
		String keyStr = result.get(condition.getRefTable().getKey());
		String[] keys = keyStr.split(",");
		for (String key : keys) {
			HashMap<String,String> annoResult = dbAdapter.getResult(condition.getRefTable().getTable(), key, condition.getFields());
			
			if (annoResult ==null || annoResult.isEmpty()){
				System.err.println("Cann't find value from table:"+condition.getRefTable().getTable()+". Key:"+key);
				return null;
			}
			
			for (String alt : condition.getAlts()) {
				results.add(alt, annoResult);
			}
		}
		
		return results;
	}

}
