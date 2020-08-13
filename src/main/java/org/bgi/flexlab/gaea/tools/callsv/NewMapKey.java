package org.bgi.flexlab.gaea.tools.callsv;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import org.apache.hadoop.io.WritableComparable;

/**
 * NewMapKey类，含有成员变量chr和pos，可作为map的输出key的类型，同时可以按照pos从小到大排序
 * @author Huifang Lu
 *
 */
public class NewMapKey implements WritableComparable<NewMapKey>{
	
	private String chr;
	private int pos;
	private int end;
	
	/**
	 * 空构造函数
	 */
	public NewMapKey() {}
	
	/**
	 * 带有两个参数（String， int）的构造函数
	 * @param chr 染色体
	 * @param pos 比对上的位置
	 */
	public NewMapKey(String chr, int pos, int end) {
		this.chr = chr;
		this.pos = pos;
		this.end = end;
	}

	public String getChr() {
		return chr;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public int getPos() {
		return pos;
	}

	public void setPos(int pos) {
		this.pos = pos;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public void readFields(DataInput in) throws IOException {
		this.chr = in.readUTF();
		this.pos = in.readInt();
		this.end = in.readInt();
	}

	public void write(DataOutput out) throws IOException {
		out.writeUTF(chr);
		out.writeInt(pos);
		out.writeInt(end);
	}
	
	@Override
	public String toString() {
		return this.chr + "\t" + this.pos + "\t" + this.end;
	}

	/**
	 * 重写了compareTo（）方法，只比较chr
	 */
	public int compareTo(NewMapKey o) {
		int num = o.getChr().compareTo(this.getChr());
		if(num==0) {
			int npos = this.pos - o.pos;
			if(npos == 0) 
				return new Integer(this.end).compareTo(new Integer(o.end));
			return npos;
		}
		return num;
	}

	@Override
	public int hashCode() {
		return chr.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(!(obj instanceof NewMapKey))
			throw new ClassCastException("Can not cast to NewMapKey class!");
		NewMapKey n = (NewMapKey) obj;
		 return this.chr.equals(n.chr) && this.pos==n.pos && this.end==n.end;
	}

	

}
