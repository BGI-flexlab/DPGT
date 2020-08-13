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

import org.bgi.flexlab.gaea.tools.annotator.AnnotationContext;
import org.bgi.flexlab.gaea.tools.annotator.VcfAnnoContext;
import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo.DbType;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;

public class DBAnnotator implements Serializable{
	
	private static final long serialVersionUID = -3944211982294335404L;
	private static HashMap<String, DBQuery> DbQueryMap = new HashMap<>();
	
	private Config config;
	private HashMap<String, Condition> dbConditionHashMap;
	
	public DBAnnotator(Config config){
		this.config = config;
		dbConditionHashMap = new HashMap<>();
	}
	
	public void annotate(VcfAnnoContext vac) throws IOException {
		
		List<String> dbNameList = config.getDbNameList();
		for (String dbName : dbNameList) {
			Condition condition = dbConditionHashMap.get(dbName);
			condition.createConditionMap(vac);
			
			DBQuery dbQuery = DbQueryMap.get(dbName);
			if (!dbQuery.executeQuery(condition)) continue;

//			alt & gene
			for (AnnotationContext annotationContext : vac.getAnnotationContexts()) {
				HashMap<String,String> result = dbQuery.getMergeResult(annotationContext);
				if (result == null) continue;
				for (String field : config.getFieldsByDB(dbName)) {
					annotationContext.putAnnoItem(field, result.get(field),false);
				}
			}
		}
	}

	public boolean annotate(VcfAnnoContext vac, String dbName) throws IOException {
		Condition condition = dbConditionHashMap.get(dbName);

		if(!dbConditionHashMap.containsKey(dbName)){
			System.err.println(dbName + " condition is null!");
			return false;
		}

		condition.createConditionMap(vac);

		DBQuery dbQuery = DbQueryMap.get(dbName);
		if (!dbQuery.executeQuery(condition)) return false;

		List<AnnotationContext> annotationContexts = new ArrayList<>();

		for(String alt: vac.getAlts()){
			LinkedList<HashMap<String,String>> resultList = dbQuery.getResultList(alt);
			for(HashMap<String,String> result: resultList){
				AnnotationContext ac = new AnnotationContext();
				ac.setGenotype(alt);
				ac.setAllele(alt);
				for (Entry<String, String> entry : result.entrySet()) {
					ac.putAnnoItem(entry.getKey(), entry.getValue(),false);
				}
				annotationContexts.add(ac);
			}
		}

		vac.setAnnotationContexts(annotationContexts);
		return true;
	}

	public boolean insert(VcfAnnoContext vac, String dbName) throws IOException {
		Condition condition = dbConditionHashMap.get(dbName);
		condition.createConditionMap(vac);
		DBQuery dbQuery = DbQueryMap.get(dbName);
		return dbQuery.insert(condition, vac.toAnnotationMaps(config.getFieldsWithoutVariant()));
	}

	public void connection() throws InstantiationException, IllegalAccessException, ClassNotFoundException, IOException {
		List<String> dbNameList = config.getDbNameList();
		List<String> dbNames = new ArrayList<>();
		dbNames.add("ANNO");
		dbNames.addAll(dbNameList);
		for (String dbName : dbNames) {
			if(!config.getDatabaseJson().hasDatabaseInfo(dbName))
				continue;
			DatabaseInfo databaseInfo = config.getDatabaseJson().getDatabaseInfo(dbName);
			String queryClassName = "DBQuery";
			if (databaseInfo.getQueryClassName() != null) {
				queryClassName = databaseInfo.getQueryClassName();
			}else {
				System.err.println("queryClassName is null, use DBQuery.class defaultly. This maybe a bug! dbName:" +dbName);
			}
			DBQuery dbQuery = (DBQuery)Class.forName("org.bgi.flexlab.gaea.tools.annotator.db." + queryClassName).newInstance();
			DbType dbType = databaseInfo.getDatabase();
			String connInfo = config.getDatabaseJson().getConnectionInfo(dbType);
			if(dbType == DbType.MYBED || dbType == DbType.VCF) {
				String tableName = config.getDatabaseJson().getDatabaseInfo(dbName).getRefTable(config.getRef()).getTable();
				dbQuery.connection(tableName, dbType, connInfo);
			}else
				dbQuery.connection(dbName, dbType, connInfo);
			DbQueryMap.put(dbName, dbQuery);

			Condition condition = new Condition(dbName,databaseInfo);

			if(dbName.equals("ANNO")){
				condition.setFields(config.getFields());
			}else
				condition.setFields(config.getFieldsByDB(dbName));
			dbConditionHashMap.put(dbName, condition);
		}
	}

	public void disconnection() throws IOException {
		Set<String> keys = DbQueryMap.keySet();
		for (String key : keys) {
			DBQuery dbQuery = DbQueryMap.get(key);
			if (dbQuery != null) {
				dbQuery.disconnection();
			}
		}
	}

}
