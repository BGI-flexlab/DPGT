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
package org.bgi.flexlab.gaea.tools.mapreduce.realigner;

import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.realigner.RealignerWriter;

import java.io.IOException;

public class RealignerContextWriter extends RealignerWriter{

	@SuppressWarnings("rawtypes")
	private Context context = null;
	private SamRecordWritable value = null;
	
	@SuppressWarnings("rawtypes")
	public RealignerContextWriter(Context context){
		this.context = context;
		this.value = new SamRecordWritable();
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public void write(GaeaSamRecord read) {
		value.set(read);
		try {
			context.write(NullWritable.get(), value);
		} catch (IOException | InterruptedException e) {
			throw new RuntimeException(e.toString());
		}
	}

	@Override
	public void close() {
	}
	
	@SuppressWarnings("rawtypes")
	public Context getContext(){
		return context;
	}
}
