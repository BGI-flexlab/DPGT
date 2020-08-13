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

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.hbase.Cell;
import org.apache.hadoop.hbase.CellUtil;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.TableName;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author huangzhibo
 *
 */
public class HbaseAdapter extends DBAdapter {

	private String columnFamily = "data";

	static Configuration conf = null;
    static Connection conn; 
    
    public HbaseAdapter(String confDir){
    	conf = HBaseConfiguration.create();
    	conf.addResource(new Path(confDir + "/hbase-site.xml"));
    	conf.addResource(new Path(confDir + "/core-site.xml"));
    }
    
    @Override
	public void connection(String cf) throws IOException {
		setColumnFamily(cf);
		conn = ConnectionFactory.createConnection(conf);
	}
	
    @Override
	public void disconnection() throws IOException {
			conn.close();
	}

    @Override
	public HashMap<String, String> getResult(String tableName, String rowKey) throws IOException {
		Get get = new Get(Bytes.toBytes(rowKey));
		get.addFamily(Bytes.toBytes(getColumnFamily()));
		Table table = conn.getTable(TableName.valueOf(tableName));
		if (!table.exists(get))
			return null;
		Result result = table.get(get);

		HashMap<String,String> resultMap = new HashMap<>();
		for  (Cell cell : result.rawCells()) {
			String key = Bytes.toString(CellUtil.cloneQualifier(cell));
			String value = Bytes.toString(CellUtil.cloneValue(cell));
			resultMap.put(key, value);
		}
		return resultMap;
	}

	public HashMap<String, String> getResult(String tableName,
			String rowKey, List<String> fields) throws IOException{
		HashMap<String,String> resultMap = new HashMap<>();
		Table table = conn.getTable(TableName.valueOf(tableName));
		Get get = new Get(Bytes.toBytes(rowKey));
		get.addFamily(Bytes.toBytes(getColumnFamily()));
		if (!table.exists(get))
			return null;
		Result result = table.get(get);
		for (String field : fields) {
			byte[] value = result.getValue(Bytes.toBytes(getColumnFamily()), Bytes.toBytes(field));
			if(value == null) continue;
			String v = Bytes.toString(value);
			if(v.contains("|"))
				v = v.replaceAll("\\|","&");
			resultMap.put(field, v);
		}
		return resultMap;
	}
	
	public boolean insert(String tableName, String rowKey, Map<String,String> fields) throws IOException{
		Table table = conn.getTable(TableName.valueOf(tableName));
		Put put = new Put(Bytes.toBytes(rowKey));
		for(String key: fields.keySet()){
			String value = fields.get(key);
			if(key == null || value == null)
				continue;
			put.addImmutable(Bytes.toBytes(getColumnFamily()), Bytes.toBytes(key), Bytes.toBytes(value));

		}
		table.put(put);
		table.close();
		return true;
	}


	public String getColumnFamily() {
		return columnFamily;
	}

	public void setColumnFamily(String columnFamily) {
		this.columnFamily = columnFamily;
	}
}
