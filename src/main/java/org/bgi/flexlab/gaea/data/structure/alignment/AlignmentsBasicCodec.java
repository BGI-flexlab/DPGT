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

import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

public class AlignmentsBasicCodec extends AlignmentsCodec<AlignmentsBasic>{

	public AlignmentsBasicCodec() {
		super();
	}
	public AlignmentsBasicCodec(Input dataInput) {
		super(dataInput);
	}

	public AlignmentsBasicCodec(Output dataOutput) {
		super(dataOutput);
	}

	@Override
	public void AlignmentsInit() {
		alignments = new AlignmentsBasic();
	}

	@Override
	protected void writeOtherInfo(AlignmentsBasic alignments) {
		// TODO Auto-generated method stub
		binaryCodec.writeInt(alignments.getSampleIndex());
		//dataOutput.flush();
	}

	@Override
	protected void readOtherInfo() {
		// TODO Auto-generated method stub
		alignments.setSampleIndex(binaryCodec.readInt());
	}
}
