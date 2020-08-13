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

import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo.DbType;

import java.io.Serializable;
import java.util.HashMap;

public class DatabaseJson implements Serializable{

	private static final long serialVersionUID = -3360762372473610852L;
	
	private final HashMap<DbType, String> connectionInfo;
	private final HashMap<String, DatabaseInfo> databaseInfo;

	public DatabaseJson(HashMap<DbType, String> connectionInfo, HashMap<String, DatabaseInfo> databaseInfo) {
		this.connectionInfo = connectionInfo;
		this.databaseInfo = databaseInfo;
	}

	public String getConnectionInfo(DbType connString) {
		return connectionInfo.get(connString);
	}

	public DatabaseInfo getDatabaseInfo(String dbName) {
		return databaseInfo.get(dbName);
	}

	public boolean hasDatabaseInfo(String dbName) {
		return databaseInfo.containsKey(dbName);
	}
}
