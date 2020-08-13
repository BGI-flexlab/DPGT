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

public abstract class GaeaFilesReader{
	protected int currentFileIndex = 0;
	protected String currentLine = null;
	
	/**
	 * get all file for input path
	 */
	public abstract void traversal(String path);
	
	/**
	 * has next line can be read?
	 */
	public abstract boolean hasNext();
	
	protected abstract int size();
	
	protected boolean filter(String name){
		if(name.startsWith("_"))
			return true;
		return false;
	}
	
	/**
	 * return next line
	 */
	public String next(){
		return currentLine;
	}
	
	/**
	 * clear all data
	 */
	public abstract void clear();
}
