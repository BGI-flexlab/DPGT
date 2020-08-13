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
package org.bgi.flexlab.gaea.util;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorDatum;

import java.util.ArrayList;
import java.util.List;

public class NestedObjectArray<T> {
	private int size;
	private int[] maximumArray;
	private Object[] data;
	
	public NestedObjectArray(final int... elements){
		size = elements.length;
		
		if(size == 0)
			throw new UserException("nested object array element size must more than 0.!");
		
		maximumArray = elements.clone();
		data = new Object[elements[0]];
	}
	
	public void put(T value,int... keys){
		if(keys.length != size)
			throw new UserException("element size must equal to "+size);
		int length = size - 1;
		
		Object[] myData = data;
		
		for(int i  = 0 ; i < length ; i++){
			if(keys[i] >= maximumArray[i])
				throw new RuntimeException("key must less than "+ maximumArray[i]);
			
			Object[] temp = (Object[]) myData[keys[i]];
			
			if(temp == null){
				temp = new Object[maximumArray[i+1]];
				myData[keys[i]] = temp;
			}
			
			myData = temp;
		}
		
		myData[keys[length]] = value;
	}
	
	@SuppressWarnings("unchecked")
	public T get(final int... keys) {
        final int length = size - 1;
        Object[] myData = data;

        for( int i = 0; i < length; i++ ) {
            if ( keys[i] >= maximumArray[i] )
                return null;
            myData = (Object[])myData[keys[i]];
            if ( myData == null )
                return null;
        }
        return (T)myData[keys[length]];
    }
	
	public List<T> getAllValues() {
        final List<T> result = new ArrayList<T>();
        fillAllValues(data, result);
        return result;
    }

    @SuppressWarnings("unchecked")
	private void fillAllValues(final Object[] array, final List<T> result) {
        for ( Object value : array ) {
            if ( value == null )
                continue;
            if ( value instanceof Object[] )
                fillAllValues((Object[])value, result);
            else
                result.add((T)value);
        }
    } 

    public static class Leave {
        public final int[] keys;
        public final Object value;

        public Leave(final int[] keys, final Object value) {
            this.keys = keys;
            this.value = value;
        }
        
        public String toString(int index){
        	StringBuilder sb = new StringBuilder();
			sb.append(index);
			sb.append("\t");
			for (int k : keys) {
				sb.append(k);
				sb.append("\t");
			}
			sb.append(((RecalibratorDatum)value).toString());
			return sb.toString();
        }
    }

    public List<Leave> getAllLeaves() {
        final List<Leave> result = new ArrayList<Leave>();
        fillAllLeaves(data, new int[0], result);
        return result;
    }

    private void fillAllLeaves(final Object[] array, final int[] path, final List<Leave> result) {
        for ( int key = 0; key < array.length; key++ ) {
            final Object value = array[key];
            if ( value == null )
                continue;
            final int[] newPath = appendToPath(path, key);
            if ( value instanceof Object[] ) {
                fillAllLeaves((Object[]) value, newPath, result);
            } else {
                result.add(new Leave(newPath, value));
            }
        }
    }

    private int[] appendToPath(final int[] path, final int newKey) {
        final int[] newPath = new int[path.length + 1];
        for ( int i = 0; i < path.length; i++ )
            newPath[i] = path[i];
        newPath[path.length] = newKey;
        return newPath;
    }
}
