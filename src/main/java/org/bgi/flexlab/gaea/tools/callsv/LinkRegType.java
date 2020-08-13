package org.bgi.flexlab.gaea.tools.callsv;

import java.util.Map;
import java.util.TreeMap;

public class LinkRegType {
	
	private String type;
	private int readNum;
	private int size;
	private Map<String, Integer> libNum;
	
	public LinkRegType() {
		this.libNum = new TreeMap<String, Integer>();
	}
	
	public LinkRegType(Reads r) {
		this.type = r.getType();
		this.libNum = new TreeMap<String, Integer>();
	}
	
	public LinkRegType(String type, int readNum, int size, Map<String, Integer> libNum) {
		this.type = type;
		this.readNum = readNum;
		this.size = size;
		this.libNum = libNum;
	}
	
	
	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public int getReadNum() {
		return readNum;
	}

	public void setReadNum(int readNum) {
		this.readNum = readNum;
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public Map<String, Integer> getLibNum() {
		return libNum;
	}

	public void setLibNum(Map<String, Integer> libNum) {
		this.libNum = libNum;
	}

	public void updateType(Reads r, int size) {
		this.readNum ++;
		this.size = this.size + size;
		
		Integer num = libNum.get(r.getLib());
		if(num == null)
			num = 0;
		num ++;
		libNum.put(r.getLib(), num);
	}

}
