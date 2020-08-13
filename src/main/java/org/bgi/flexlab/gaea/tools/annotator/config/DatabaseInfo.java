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

import java.io.Serializable;
import java.util.HashMap;

public class DatabaseInfo implements Serializable {

	public enum DbType {
		HBASE
		, MYSQL
		, MYBED
		,VCF
		,TSV
	}
	
	public static enum ConditionKey {
		CHR
		, CHROM
		, START  //0-base, POS - 1
		, POS
		, END0   //0-base
		, END  //1-base , some database use this , example: HGMD
		, REF
		, ALT
		, GENE
		, ASSEMBLY
	}
	
	private static final long serialVersionUID = -5850759043237677551L;
	
	private final DbType database;
	private final String queryClassName;
	private final String queryCondition;
	private final String altField;
	
	private HashMap<String, String> fields;
	private HashMap<String, RefTableInfo> dbVersion;
	
	public DatabaseInfo(DbType database, String queryCondition, String queryClassName, String extendField,
			String altField, HashMap<String, RefTableInfo> dbVersion, HashMap<String, String> fields) {
		this.database = database;
		this.queryClassName = queryClassName;
		this.queryCondition = queryCondition;
		this.altField = altField;
		this.fields = fields;
		this.dbVersion = dbVersion;
	}
	
	public void setFields(HashMap<String, String> fields) {
		this.fields = fields;
	}

	public HashMap<String, String> getFields() {
		if (fields == null) {
			fields = new HashMap<>();
		}
		return fields;
	}

	public String getAltField() {
		return altField;
	}

	/**
	 * key is reference version for human : GRCh37, GRCh38
	 * @param key
	 * @return
	 */
	public RefTableInfo getRefTable(String key) {
		return dbVersion.get(key);
	}

	public String getQueryCondition() {
		return queryCondition;
	}

	public String getQueryClassName() {
		return queryClassName;
	}

	public DbType getDatabase() {
		return database;
	}

}
