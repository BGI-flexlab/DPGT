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
package org.bgi.flexlab.gaea.util;

public class ChromosomeUtils {
	
	public static String formatChrName(String chrName) {
		chrName = chrName.toLowerCase();
		if(!chrName.startsWith("chr")) {
			chrName = "chr" + chrName;
		}
		return chrName;
	}

	public static String getNoChrName(String chrName) {
		if(chrName.startsWith("chr"))
			return chrName.substring(3);
		else
			return chrName;
	}
}
