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
import org.bgi.flexlab.gaea.tools.annotator.config.Config;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo.ConditionKey;
import org.bgi.flexlab.gaea.tools.annotator.config.DatabaseInfo.DbType;
import org.bgi.flexlab.gaea.tools.annotator.config.RefTableInfo;

import java.io.Serializable;
import java.util.*;

/**
 * @author huangzhibo
 *
 * 获取不同数据库的查询条件
 */
public class Condition implements Serializable{

	private static final long serialVersionUID = 5474372956720082769L;
	
	private String conditionString = null;
	private String dbName;
	private DbType dbType;   //数据库类型: Hbase, Mysql
	private List<String> fields;
	private List<String> alts;
	private Set<String> genes;
	private Set<String> transcriptIds;
	private DatabaseInfo dbInfo;
	private RefTableInfo refTable;
	
	private boolean keyHasGene;
	private boolean keyHasAlt;
	private boolean keyIsChrPos;
	
	
	private LinkedHashMap<ConditionKey, String> conditionMap;
	
	public Condition(String dbName, DatabaseInfo dbInfo){
		this.dbInfo = dbInfo;
		this.dbType = dbInfo.getDatabase();
		setKeyIsChrPos();
		conditionMap = new LinkedHashMap<ConditionKey, String>();
		this.refTable = dbInfo.getRefTable(Config.getConfigInstance().getRef());
	}
	
	public Condition(DatabaseInfo dbInfo, LinkedHashMap<ConditionKey, String>conditionMap){
		this.conditionMap = conditionMap;
		this.dbInfo = dbInfo;
		this.setRefTable(dbInfo.getRefTable(conditionMap.get(ConditionKey.ASSEMBLY)));
		this.dbType = dbInfo.getDatabase();
	}
	
	public Condition(String dbName, DatabaseInfo dbInfo, LinkedHashMap<ConditionKey, String>conditionMap){
		setDbName(dbName);
		this.conditionMap = conditionMap;
		this.dbInfo = dbInfo;
		this.setRefTable(dbInfo.getRefTable(conditionMap.get(ConditionKey.ASSEMBLY)));
		this.dbType = dbInfo.getDatabase();
	}

	/**
	 * condition string type
	 */
	public String getConditionString() {
		StringBuilder sb = new StringBuilder();
		String[] conKeys = dbInfo.getQueryCondition().split("_");
		if(dbType == DbType.HBASE){
			for (int i = 0; i < conKeys.length-1; i++) {
				sb.append(conditionMap.get(ConditionKey.valueOf(conKeys[i]))+'-');
			}
			sb.append(conditionMap.get(ConditionKey.valueOf(conKeys[conKeys.length-1])));
		}else if(dbType == DbType.VCF){
			for (int i = 0; i < conKeys.length-1; i++) {
				sb.append(conditionMap.get(ConditionKey.valueOf(conKeys[i]))+'\t');
			}
			sb.append(conditionMap.get(ConditionKey.valueOf(conKeys[conKeys.length-1])));
		}else if(dbType == DbType.MYSQL){
			sb.append("where ");
			for (int i = 0; i < conKeys.length-1; i++) {
				sb.append(conKeys[i] + " = " + conditionMap.get(conKeys[i])+" and ");
			}
			sb.append(conKeys[conKeys.length-1] + " = " + conditionMap.get(conKeys[conKeys.length-1]));
		}else {
			for (int i = 0; i < conKeys.length-1; i++) {
				sb.append(conditionMap.get(ConditionKey.valueOf(conKeys[i]))+'-');
			}
			sb.append(conditionMap.get(ConditionKey.valueOf(conKeys[conKeys.length-1])));
		}
		conditionString = sb.toString();
		return conditionString;
	}
	
	/**
	 * condition string hash
	 */
	public HashMap<String,String> getConditionHash() {
		HashMap<String,String> conditionHash = new HashMap<String, String>();
		if(dbInfo.getQueryCondition().indexOf("GENE") != -1){
			for (String gene: genes) {
				conditionMap.put(ConditionKey.GENE, gene);
				String conditionStr = getConditionString();
				conditionHash.put(gene, conditionStr);
			}
		}else if (dbInfo.getQueryCondition().indexOf("ALT") != -1) {
			for (String alt: alts) {
				conditionMap.put(ConditionKey.ALT, alt);
				String conditionStr = getConditionString();
				conditionHash.put(alt, conditionStr);
			}
		}else {
			conditionHash.put("common", getConditionString());
		}
		
		return conditionHash;
	}
	
	/**
	 * condition string list
	 */
	@Deprecated
	public List<String> getConditionList() {
		List<String> conditionList = new ArrayList<String>();
		if(dbInfo.getQueryCondition().indexOf("GENE") != -1){
			for (String gene: genes) {
				conditionMap.put(ConditionKey.GENE, gene);
				String conditionStr = getConditionString();
				conditionList.add(conditionStr);
			}
		}else if (dbInfo.getQueryCondition().indexOf("ALT") != -1) {
			for (String alt: alts) {
				conditionMap.put(ConditionKey.ALT, alt);
				String conditionStr = getConditionString();
				conditionList.add(conditionStr);
			}
		}else {
			conditionList.add(getConditionString());
		}
		
		return conditionList;
	}
	
	public List<String> getConditionAltList() {
		List<String> conditionList = new ArrayList<String>();
		
		for (String alt: alts) {
			conditionMap.put(ConditionKey.ALT, alt);
			String conditionStr = getConditionString();
			conditionList.add(conditionStr);
		}
		
		return conditionList;
	}
	
	public String getDbName() {
		return dbName;
	}

	public void setDbName(String dbName) {
		this.dbName = dbName;
	}
	
	public List<String> getAlts() {
		return alts;
	}

	public void setAlts(List<String> alts) {
		this.alts = alts;
	}

	public DbType getDbType() {
		return dbType;
	}

	public void setDbType(DbType dbType) {
		this.dbType = dbType;
	}

	public DatabaseInfo getTableInfo() {
		return dbInfo;
	}

	public void setDbInfo(DatabaseInfo dbInfo) {
		this.dbInfo = dbInfo;
	}

	public void setFields(List<String> fields) {
		this.fields = fields;
	}

	public List<String> getFields(){
		return fields;
	}

	public String getAltField(){
		return dbInfo.getAltField();
	}

	public void setConditionMap(LinkedHashMap<ConditionKey, String> conditionMap) {
		this.conditionMap = conditionMap;
	}

	public HashMap<ConditionKey, String> getConditionMap() {
		return conditionMap;
	}

	public RefTableInfo getRefTable() {
		return refTable;
	}

	public void setRefTable(RefTableInfo refTable) {
		this.refTable = refTable;
	}

	public void setGenes(Set<String> genes) {
		this.genes = genes;
	}

	public Set<String> getGenes() {
		return genes;
	}

	public Set<String> getTranscriptIds() {
		return transcriptIds;
	}

	public void setTranscriptIds(Set<String> transcriptIds) {
		this.transcriptIds = transcriptIds;
	}

	public boolean isKeyHasGene() {
		return keyHasGene;
	}

	public void setKeyHasGene(boolean keyHasGene) {
		this.keyHasGene = keyHasGene;
	}

	public boolean isKeyHasAlt() {
		return keyHasAlt;
	}

	public void setKeyHasAlt(boolean keyHasAlt) {
		this.keyHasAlt = keyHasAlt;
	}

	public boolean isKeyIsChrPos() {
		return keyIsChrPos;
	}

	private void setKeyIsChrPos() {
		if (!keyHasAlt && !keyHasGene) {
			this.keyIsChrPos = true;
		}else {
			this.keyIsChrPos = false;
		}
	}

	public void createConditionMap(VcfAnnoContext vac) {
		// TODO Auto-generated method stub
		conditionMap.put(ConditionKey.CHR, vac.getChromeNoChr());
		conditionMap.put(ConditionKey.CHROM, vac.getChrome());
		conditionMap.put(ConditionKey.ASSEMBLY, Config.getConfigInstance().getRef());
		conditionMap.put(ConditionKey.POS, String.valueOf(vac.getStart()));
		conditionMap.put(ConditionKey.START, String.valueOf(vac.getStart()-1));
		conditionMap.put(ConditionKey.END, String.valueOf(vac.getEnd()));
		setAlts(vac.getAlts());  // alt == annotationContext.getAllele() ?? TODO
		setTranscriptIds(vac.getTranscriptIds());
		setGenes(vac.getGenes());
	}
	
}
