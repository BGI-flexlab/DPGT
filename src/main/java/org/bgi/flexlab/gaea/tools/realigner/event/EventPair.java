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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.realigner.event;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

import java.util.TreeSet;

public class EventPair {
	private Event left,right;
	private TreeSet<GenomeLocation> intervals = new TreeSet<GenomeLocation>();
	
	public EventPair(Event l,Event r){
		this.left = l;
		this.right = r;
	}
	
	public EventPair(Event l,Event r,TreeSet<GenomeLocation> set1,TreeSet<GenomeLocation> set2){
		this.left = l;
		this.right = r;
		intervals.addAll(set1);
		intervals.addAll(set2);
	}
	
	public TreeSet<GenomeLocation> getIntervals(){
		return intervals;
	}
	
	public Event getLeft(){
		return left;
	}
	
	public Event getRight(){
		return right;
	}
	
	public void setLeft(Event l){
		left = l;
	}
	
	public void setRight(Event r){
		right = r;
	}
	
	public String toString(){
    	StringBuilder sb = new StringBuilder();
    	
    	if(left != null && right == null)
    		sb.append(left.toString()+";");
    	else if(left != null && right != null){
    		sb.append(left.toString()+";"+right.toString());
    	}
    	return sb.toString();
    }
}
