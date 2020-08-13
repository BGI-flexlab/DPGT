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
package org.bgi.flexlab.gaea.tools.mapreduce.annotator;


import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.data.mapreduce.writable.PairWritable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class AnnoSortMapper extends Mapper<LongWritable, Text, PairWritable, Text> {

	private PairWritable resultKey;
	private Text resultValue;

	@Override
	protected void setup(Context context)
			throws IOException, InterruptedException {
		resultKey = new PairWritable();
		resultValue = new Text();
	}

	@Override
	protected void map(LongWritable key, Text value, Context context)
			throws IOException, InterruptedException {
		String annoLine = value.toString();
		if (annoLine.startsWith("#")) return;

		String[] fields = annoLine.split("\t", 4);
		String secondKey = fields[1] + "-" + String.format("%09d",Integer.parseInt(fields[2]));

		resultKey.set(fields[0], secondKey);
		List<String> vcfField = new ArrayList<>();
		vcfField.add(fields[1]);
		vcfField.add(fields[2]);
		vcfField.add(fields[3]);
		resultValue.set(String.join("\t",vcfField));
		context.write(resultKey, resultValue);
	}
	
	@Override
	protected void cleanup(Context context)
			throws IOException, InterruptedException {

	}
}
