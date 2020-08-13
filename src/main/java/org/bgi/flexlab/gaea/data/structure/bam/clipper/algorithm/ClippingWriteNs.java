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
package org.bgi.flexlab.gaea.data.structure.bam.clipper.algorithm;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.bam.clipper.ClippingRegion;

public class ClippingWriteNs implements ReadClippingAlgorithm{

	@Override
	public GaeaSamRecord apply(GaeaSamRecord read, ClippingRegion cr) {
		byte[] bases = read.getReadBases();
		
		for(int i = cr.getStart() ; i < cr.getStop() ; i++){
			bases[i] = 'N';
		}
		
		read.setReadBases(bases);
		
		return read;
	}
}
