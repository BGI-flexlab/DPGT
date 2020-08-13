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

import java.util.*;
import java.util.Map.Entry;

public class DbSNPQuery extends DBQuery {

	@Override
	protected List<HashMap<String, String>> splitResult(HashMap<String, String> result, int altNum) {
		List<HashMap<String, String>> resultList = new ArrayList<>();
		if(altNum == 1) {
			if(result.containsKey("CAF")){
				String cafValue = result.get("CAF");
				String[] caf = cafValue.split(",");
				if(caf.length == 2)
					result.put("CAF", caf[1]);
			}
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
				if(altNum + 1 == values.length){
					for (int i = 0; i < altNum; i++) {
						resultList.get(i).put(entry.getKey(), values[i+1]);
					}
				}else if(altNum == values.length){
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
}
