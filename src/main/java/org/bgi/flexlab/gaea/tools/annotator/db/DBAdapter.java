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
import java.util.List;
import java.util.Map;

public abstract class DBAdapter {
	
	public abstract void connection(String dbName) throws IOException;
	
	public abstract void disconnection() throws IOException;

	/**
	 * 根据查询条件和fieldMap查询所有相关字段，fieldMap的value为数据库中字段String，其key值作为返回的HashMap的key。
	 *
	 * @param tableName
	 * @param condition
	 * @param fields
	 * @return
	 * @throws IOException
	 */
	public HashMap<String, String> getResult(String tableName, String condition, List<String> fields) throws IOException {
		return null;
	}

	public List<HashMap<String, String>> getResult(String condition, List<String> fields) throws IOException {
		return null;
	}

	/**
	 * 根据查询条件查询所有相关信息，返回的HashMap包含该行所有字段
	 *
	 * @return HashMap<String, String>
	 * @tableName 表名
	 * @conditionString conditionString (对于Hbase则为rowKey)
	 */
	public HashMap<String, String> getResult(String tableName,
											 String conditionString) throws IOException {
		return null;
	}

	public boolean insert(String tableName, String rowKey, Map<String, String> fields) throws IOException {
		return false;
	}

}
