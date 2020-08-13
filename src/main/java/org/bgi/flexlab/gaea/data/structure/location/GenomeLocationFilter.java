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
package org.bgi.flexlab.gaea.data.structure.location;

import org.bgi.flexlab.gaea.util.Window;

import java.util.ArrayList;

public abstract class GenomeLocationFilter {
	public abstract boolean filter(GenomeLocation location,Window win);
	
	public ArrayList<GenomeLocation> filterList(ArrayList<GenomeLocation> intervals,Window win){
		ArrayList<GenomeLocation> filtered = new ArrayList<GenomeLocation>();
		
		for(GenomeLocation location : intervals){
			if(!filter(location,win)){
				filtered.add(location);
			}
		}
		
		return filtered;
	}
}
