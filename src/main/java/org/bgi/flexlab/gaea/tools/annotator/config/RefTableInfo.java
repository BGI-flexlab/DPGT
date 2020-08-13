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
package org.bgi.flexlab.gaea.tools.annotator.config;

public class RefTableInfo {
	private final String table;
	private final String indexTable;
	private final String key;
	
	public RefTableInfo(String table, String indexTable, String key) {
		this.table = table;
		this.indexTable = indexTable;
		this.key = key;
	}

	public String getTable() {
		return table;
	}

	public String getIndexTable() {
		return indexTable;
	}

	public String getKey() {
		return key;
	}

}
