package org.bgi.flexlab.gaea.tools.callsv;

import java.util.ArrayList;
import java.util.List;


public class Reads {

	private String name;
	private String chr;
	private String type;
	private String lib;
	private int insert;
	private List<Integer> reg;
	
	public Reads(){
		setDefault();
	}
	
	public Reads(String name, String chr, String type, String lib, int insert, List<Integer> reg){
		this.name = name;
		this.chr = chr;
		this.type = type;
		this.lib = lib;
		this.insert = insert;
		this.reg = reg;
	}
	
	public Reads(SamWritable r) {
		this.name = r.getReadName();
		this.chr = r.getChr();
		this.type = r.getType();
		this.lib = r.getLib();
		this.insert = r.getInsert();
		this.reg = new ArrayList<Integer>();
	}
	
	private void setDefault() {
		this.name = null;
		this.chr = null;
		this.type = null;
		this.lib = null;
		this.insert = 0;
		this.reg = new ArrayList<Integer>();
	}
	
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getChr() {
		return chr;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public String getType() {
		return type;
	}
	
	public void setType(String type) {
		this.type = type;
	}
	
	public String getLib() {
		return lib;
	}
	
	public void setLib(String lib) {
		this.lib = lib;
	}
	
	public int getInsert() {
		return insert;
	}

	public void setInsert(int insert) {
		this.insert = insert;
	}

	public List<Integer> getReg() {
		return reg;
	}
	
	public void setReg(List<Integer> reg) {
		this.reg = reg;
	}

	@Override
	public int hashCode() {
		return this.name.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(!(obj instanceof Reads))
			throw new ClassCastException("obj can not cast to Reads class!");
		Reads r = (Reads)obj;
		return this.name.equals(r.name);
	}

	@Override
	public String toString() {
		return name + "\t" + chr + "\t" + type + "\t" + lib + "\t" + insert +"\t"+ reg ;
	}
	
}
