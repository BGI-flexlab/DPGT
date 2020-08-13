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
package org.bgi.flexlab.gaea.data.structure.alignment;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.RuntimeEOFException;
import org.bgi.flexlab.gaea.data.structure.bam.SAMCompressionInformationBasic;

import java.io.InputStream;
import java.io.OutputStream;

public abstract class AlignmentsCodec <T extends SAMCompressionInformationBasic>{
	protected BinaryCodec binaryCodec;
	protected T alignments;
	//FieldAccess access = FieldAccess.get(AlignmentsBasic.class);

	public AlignmentsCodec() {
		binaryCodec = new BinaryCodec();
	}

	public AlignmentsCodec(InputStream dataInput) {
		binaryCodec = new BinaryCodec();
		this.binaryCodec.setInputStream(dataInput);
	}
	
	public AlignmentsCodec(OutputStream dataOutput) {
		binaryCodec = new BinaryCodec();
		this.binaryCodec.setOutputStream(dataOutput);
		//AlignmentsInit();
	}

	public void setInputStream(InputStream dataInput) {
		this.binaryCodec.setInputStream(dataInput);
	}

	public void setOutputStream(OutputStream dataOutput) {
		this.binaryCodec.setOutputStream(dataOutput);
	}

	public void encode(T alignments) {
		writeBasic(alignments);
		writeOtherInfo(alignments);
	}

	public T 	decode() {
		AlignmentsInit();
		if(readBasic()) {
			readOtherInfo();
			return alignments;
		} else
			return null;
	}
	
	private void writeBasic(T alignments) {
		binaryCodec.writeShort((short)alignments.getFlag());
		binaryCodec.writeInt(alignments.getChrNameIndex());
		binaryCodec.writeInt(alignments.getPosition());
		binaryCodec.writeShort(alignments.getMappingQual());
		binaryCodec.writeInt(alignments.getCigarsLength());
		int[] cigars = alignments.getCigars();
		for(int i = 0; i < alignments.getCigarsLength(); i++) {
			binaryCodec.writeInt(cigars[i]);
		}
		binaryCodec.writeInt(alignments.getReadLength());
		binaryCodec.writeBytes(alignments.getQualities());
		binaryCodec.writeBytes(alignments.getCompressedReadBasesBytes());
		//System.out.println("write record : " + alignments.toString());
	}
	
	private boolean readBasic() {
		//System.out.print("read record : ");
		try {
			alignments.setFlag(binaryCodec.readShort());
		} catch (RuntimeEOFException e) {
			return false;
		}
		alignments.setChrNameIndex(binaryCodec.readInt());
		alignments.setPosition(binaryCodec.readInt());
		alignments.setMappingQual(binaryCodec.readShort());
		int cigarLength = binaryCodec.readInt();
		int[] cigars = new int[cigarLength];
		for(int i = 0; i < cigarLength; i++) {
			cigars[i] = binaryCodec.readInt();
		}
		alignments.setCigars(cigars);
		int readLength = binaryCodec.readInt();
		int baseLength = (readLength + 1) / 2;
		byte[] qualities = new byte[readLength];
		byte[] readBases = new byte[baseLength];
		binaryCodec.readBytes(qualities);
		binaryCodec.readBytes(readBases);
		alignments.setQualities(qualities);
		alignments.getReadBasesFromCompressedReadBasesBytes(readBases, readLength);
		//System.out.println(alignments.toString());
		return true;
	}
	
	public abstract void AlignmentsInit();

	protected abstract void writeOtherInfo(T alignments);
	
	protected abstract void readOtherInfo();
	
	
}
