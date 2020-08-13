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
package org.bgi.flexlab.gaea.tools.realigner.alternateconsensus;

import htsjdk.samtools.Cigar;

import java.util.ArrayList;
import java.util.Arrays;

import org.bgi.flexlab.gaea.util.Pair;

public class AlternateConsensus {
	private final byte[] str;
	private final ArrayList<Pair<Integer, Integer>> readIndexes;
	private final int positionOnReference;
	private int mismatchSum;
	private Cigar cigar;

    public AlternateConsensus(byte[] str, Cigar cigar, int positionOnReference) {
        this.str = str;
        this.cigar = cigar;
        this.positionOnReference = positionOnReference;
        mismatchSum = 0;
        readIndexes = new ArrayList<Pair<Integer, Integer>>();
    }

    @Override
    public boolean equals(Object o) {
    	if(o instanceof AlternateConsensus){
    		AlternateConsensus other = (AlternateConsensus)o;
    		return ( this == other || Arrays.equals(this.str,other.getSequence()) ) ;
    	}
    	return false;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(this.str);
    }
    
    public byte[] getSequence(){
    	return this.str;
    }
    
    public ArrayList<Pair<Integer, Integer>> getReadIndexes(){
    	return this.readIndexes;
    }
    
    public void add(Pair<Integer,Integer> pair){
    	this.readIndexes.add(pair);
    }
    
    public int getPositionOnReference(){
    	return this.positionOnReference;
    }
    
    public int getMismatch(){
    	return this.mismatchSum;
    }
    
    public void addMismatch(int score){
    	this.mismatchSum += score;
    }
    
    public Cigar getCigar(){
    	return this.cigar;
    }
    
    public void setCigar(Cigar cigar){
    	this.cigar = cigar;
    }
    
    public void clear(){
    	this.readIndexes.clear();
    }
}
