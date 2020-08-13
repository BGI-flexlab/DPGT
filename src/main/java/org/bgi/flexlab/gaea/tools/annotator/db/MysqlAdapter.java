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
import java.sql.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class MysqlAdapter extends DBAdapter {
	
	
	private static final String USER="root";
	private static final String PASSWORD="";
	
	private static String driver ="com.mysql.jdbc.Driver";
	private static String url=null;
	private static Connection conn=null;
	
	MysqlAdapter(String url){
		MysqlAdapter.url = url;
	}
	
	@Override
	public void connection(String dbName) throws IOException {
		try {
			//1.加载驱动程序
			Class.forName(driver);
			//2.获得数据库的连接
			setConn(DriverManager.getConnection(url+dbName, USER, PASSWORD));
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
	}

	private void setConn(Connection connection) {
		conn = connection;
	}
	
	public static Connection getConn() {
		return conn;
	}

	@Override
	public void disconnection() {
		try {
			conn.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	public static String getDriver() {
		return driver;
	}
	
	public static void setDriver(String driver) {
		MysqlAdapter.driver = driver;
	}
	
	public HashMap<String, String> getResult(String tableName, String condition, List<String> tags) {
		HashMap<String,String> resultMap = new HashMap<String,String>();
		StringBuilder sb=new StringBuilder();
		sb.append("select * from " + tableName + " ");
		sb.append(condition);
		PreparedStatement ptmt = null;
		ResultSet rs = null;
		try {
			ptmt=conn.prepareStatement(sb.toString());
			rs=ptmt.executeQuery();
			for (String key : tags) {
					resultMap.put(key, rs.getString(key));
			}
		} catch (SQLException e) {
			System.err.println("SQLException: query tableName(" + tableName + ") is Wrong!\nStatement: "+ sb.toString());
			e.printStackTrace();
		}
		return resultMap;
	}
	
	@Override
	public HashMap<String, String> getResult(String indexTable,
			String conditionString) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean insert(String tableName, String rowKey, Map<String, String> fields) throws IOException {
		return false;
	}


}
