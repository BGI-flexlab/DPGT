package org.bgi.flexlab.gaea.tools.callsv;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import org.apache.hadoop.io.Writable;

import htsjdk.samtools.SAMRecord;

public class SamWritable implements Writable{
	
	private String lib;
	private String readName;
	private int flag;
	private String chr;
	private	int start;
	private int end;
	private	int insert;
	private String strand;
	private String type;
	private int readLen;
	
	/**
	 * 空构造函数
	 */
	public SamWritable(){}
	
	/**
	 * 带有String类型参的数构造函数
	 * @param line String类，包含了必要的比对信息，用制表符分开，并且顺序一致
	 */
	public SamWritable(String line){
		String[] tmp = line.split("\\s+");
		this.lib = tmp[0];
		this.readName = tmp[1];
		this.flag = Integer.parseInt(tmp[2]);
		this.chr = tmp[3];
		this.start = Integer.parseInt(tmp[4]);
		this.end = Integer.parseInt(tmp[5]);
		this.insert = Math.abs(Integer.parseInt(tmp[6]));
		this.strand = tmp[7];
		this.type = tmp[8];
		this.readLen = Integer.parseInt(tmp[9]);
	}
	
	
	/**
	 * 带有两个参数（SAMRecord类型和String类型）的构造函数
	 * @param r SAMRecord类，包含reads完整的比对信息
	 * @param type 可能的APRs类型
	 */
	public SamWritable(SAMRecord r, String type) {
		this.lib = r.getReadGroup().getLibrary();
		this.readName = r.getReadName();
		this.flag = r.getFlags();
		this.chr = r.getReferenceName();
		this.start = r.getAlignmentStart();
		this.end = r.getAlignmentEnd();
		this.insert = Math.abs(r.getInferredInsertSize());
		this.strand = (r.getReadNegativeStrandFlag() == true) ? "-" : "+";
		this.type = type;
		this.readLen = r.getReadLength();
		
	}
	
	/**
	 * 带有多个参数的构造函数
	 * @param lib String类型，文库名
	 * @param readName String类型，read名字
	 * @param flag int类型，比对的flag
	 * @param chr String类型，比对上的染色体号
	 * @param start int类型，比对上的起始位置
	 * @param end int类型，比对上的终止位置
	 * @param insert int类型，一对reads比对上的距离
	 * @param strand String类型，比对上的方向
	 * @param type String类型，APRs的类型
	 * @param readLen int类型，read的长度
	 */
	
	public SamWritable(String lib, String readName, int flag, String chr, int start, int end, int insert, String strand, String type, int readLen){
		this.lib = lib;
		this.readName = readName;
		this.flag = flag;
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.insert = Math.abs(insert);
		this.strand = strand;
		this.type = type;
		this.readLen = readLen;
	}
	
	/**
	 * set方法，带有两个参数，根据参数设定各个成员变量的值
	 * @param r SAMRecord类，包含reads完整的比对信息
	 * @param type 可能的APRs类型
	 */
	
	public void set(SAMRecord r, String type) {
		this.lib = r.getReadGroup().getLibrary();
		this.readName = r.getReadName();
		this.flag = r.getFlags();
		this.chr = r.getReferenceName();
		this.start = r.getAlignmentStart();
		this.end = r.getAlignmentEnd();
		this.insert = Math.abs(r.getInferredInsertSize());
		this.strand = (r.getReadNegativeStrandFlag() == true) ? "-" : "+";
		this.type = type;
		this.readLen = r.getReadLength();
	}
	
	/**
	 * 带有多个参数的set方法，根据参数设定或者修改各个成员变量的值
	 * @param lib String类型，文库名
	 * @param readName String类型，read名字
	 * @param flag int类型，比对的flag
	 * @param chr String类型，比对上的染色体号
	 * @param start int类型，比对上的起始位置
	 * @param end int类型，比对上的终止位置
	 * @param insert int类型，一对reads比对上的距离
	 * @param strand String类型，比对上的方向
	 * @param type String类型，APRs的类型
	 * @param readLen int类型，read的长度
	 */
	
	public void set(String lib, String readName, int flag, String chr, int start, int end, int insert, String strand, String type, int readLen) {
		this.lib = lib;
		this.readName = readName;
		this.flag = flag;
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.insert = Math.abs(insert);
		this.strand = strand;
		this.type = type;
		this.readLen = readLen;
	}
	
	/**
	 * 带有一个String类型参数的set方法
	 * @param line String类，包含了必要的比对信息，用制表符分开，并且顺序一致
	 */
	
	public void set(String line){
		String[] tmp = line.split("\\s+");
		this.lib = tmp[0];
		this.readName = tmp[1];
		this.flag = Integer.parseInt(tmp[2]);
		this.chr = tmp[3];
		this.start = Integer.parseInt(tmp[4]);
		this.end = Integer.parseInt(tmp[5]);
		this.insert = Math.abs(Integer.parseInt(tmp[6]));
		this.strand = tmp[7];
		this.type = tmp[8];
		this.readLen = Integer.parseInt(tmp[9]);
	}
	

	/**
	 * 获取文库名
	 * @return lib 文库名
	 */
 	public String getLib() {
		return lib;
	}
 	
 	/**
 	 * 设定文库名
 	 * @param lib 文库名
 	 */

	public void setLib(String lib) {
		this.lib = lib;
	}

	public String getReadName() {
		return readName;
	}

	public void setReadName(String readName) {
		this.readName = readName;
	}

	public int getFlag() {
		return flag;
	}

	public void setFlag(int flag) {
		this.flag = flag;
	}

	public String getChr() {
		return chr;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public int getInsert() {
		return insert;
	}

	public void setInsert(int insert) {
		this.insert = insert;
	}

	public String getStrand() {
		return strand;
	}

	public void setStrand(String strand) {
		this.strand = strand;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public int getReadLen() {
		return readLen;
	}

	public void setReadLen(int readLen) {
		this.readLen = readLen;
	}

	public void readFields(DataInput in) throws IOException {
		this.lib = in.readUTF();
		this.readName = in.readUTF();
		this.flag = in.readInt();
		this.chr = in.readUTF();
		this.start = in.readInt();
		this.end = in.readInt();
		this.insert = in.readInt();
		this.strand = in.readUTF();
		this.type = in.readUTF();
		this.readLen = in.readInt();
		
	}
	

	public void write(DataOutput out) throws IOException {
		out.writeUTF(lib);
		out.writeUTF(readName);
		out.writeInt(flag);
		out.writeUTF(chr);
		out.writeInt(start);
		out.writeInt(end);
		out.writeInt(insert);
		out.writeUTF(strand);
		out.writeUTF(type);
		out.writeInt(readLen);
		
	}

	
	@Override
	public String toString() {
		return lib + "\t" + readName + "\t" + flag + "\t" + chr + "\t" + start + "\t" + end + "\t" + insert + "\t" + strand + "\t" + type + "\t" + readLen;
	}


}
