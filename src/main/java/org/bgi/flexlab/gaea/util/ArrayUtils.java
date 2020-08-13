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

public class ArrayUtils {
	public static String[] subArray(String[] args,int start,int length){
		String[] sub = new String[length];
		for(int i = 0; i < length ; i++){
			sub[i] = args[start+i];
		}
		return sub;
	}
	
	public static String[] subArray(String[] args,int start){
		return subArray(args,start,args.length-start);
	}
}
