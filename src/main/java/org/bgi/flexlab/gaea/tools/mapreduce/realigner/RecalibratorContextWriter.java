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
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer.Context;
import org.apache.hadoop.mapreduce.lib.output.MultipleOutputs;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.tools.realigner.RealignerWriter;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorWriter;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTable;

import java.io.IOException;

public class RecalibratorContextWriter extends RealignerWriter implements RecalibratorWriter {
	public final static String RECALIBRATOR_TABLE_TAG = "bqsr";
	
	@SuppressWarnings("rawtypes")
	private Context context = null;
	private MultipleOutputs<NullWritable, Text> mos = null;
	private SamRecordWritable value = null;

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public RecalibratorContextWriter(Context ctx,boolean multiple) {
		if(multiple)
			mos = new MultipleOutputs<NullWritable, Text>(ctx);
		this.context = ctx;
		value = new SamRecordWritable();
	}
	
	public RecalibratorContextWriter(@SuppressWarnings("rawtypes") Context context){
		this(context,false);
	}

	@SuppressWarnings("unchecked")
	@Override
	public void write(RecalibratorTable table) {
		for (String str : table.valueStrings()) {
			try {
				if(mos == null){
					context.write(NullWritable.get(), new Text(str));
				}else{
					mos.write(RECALIBRATOR_TABLE_TAG, NullWritable.get(), new Text(str));
				}
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			} catch (InterruptedException e) {
				throw new RuntimeException(e.toString());
			}
		}
	}

	@Override
	public void close() {
		if(mos != null){
			try {
				mos.close();
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			} catch (InterruptedException e) {
				throw new RuntimeException(e.toString());
			}
		}
	}

	@SuppressWarnings("unchecked")
	@Override
	public void write(GaeaSamRecord read) {
		value.set(read);
		try {
			if(context != null)
				context.write(NullWritable.get(), value);
		} catch (IOException | InterruptedException e) {
			throw new RuntimeException(e.toString());
		}
	}
	
	@SuppressWarnings("rawtypes")
	public Context getContext(){
		return this.context;
	}
}
