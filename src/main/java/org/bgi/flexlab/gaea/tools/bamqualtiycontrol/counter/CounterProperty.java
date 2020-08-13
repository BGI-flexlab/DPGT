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
package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter;

public interface CounterProperty {
	
	public enum Interval implements CounterProperty{
		WHOLEGENOME, TARGET, FLANK, CHRX, CHRY;
	}
	
	public enum Depth implements CounterProperty{
		ABOVE_ZREO(0), FOURX(4), TENX(10), TWENTYX(20), THIRTYX(30), FIFTYX(50), HUNDREDX(100),
		TOTALDEPTH(-1);
		
		private int depth;
		
		Depth(int depth) {
			this.depth = depth;
		}
		
		public int getDepth() {
			return depth;
		}
	}
	
	public enum DepthType implements CounterProperty{
		WITHOUT_PCR(0), NORMAL(0);
		
		private int depth;
		
		private DepthType(int depth) {
			// TODO Auto-generated constructor stub
			this.depth = depth;
		}
		
		public DepthType setDepth(int depth) {
			this.depth = depth;
			return this;
		}
		
		public int getDepth() {
			return depth;
		}
	}
	
	public enum BaseType implements CounterProperty{
		TOTALBASE(1), MISMATCH(1), FOURXCOVERED(1), FOURXNONNCOVERED(1), NONNCOVERED(1), COVERED(1), TENXCOVERED(1), TENXNONNCOVERED(1), THIRTYXCOVERED(1), THIRTYXNONNCOVERED(1), INDELREF(1), MISMATCHREF(1);

		private int count;

		private BaseType(int count) {
			this.count = count;
		}
		
		public BaseType setCount(int count) {
			this.count = count;
			return this;
		}
		
		public int getCount() {
			return count;
		}
	}
	
	public enum ReadType implements CounterProperty{
		TOTALREADS, MAPPED, UNIQUE, DUP, PE, CLIPPED, MISMATCHREADS, INDEL;
	}
}
