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
package org.bgi.flexlab.gaea.tools.recalibrator.table;

import htsjdk.samtools.SAMFileHeader;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFilesReader;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorContextWriter;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorDatum;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.CovariateUtil;
import org.bgi.flexlab.gaea.util.NestedObjectArray;

public class RecalibratorTableCombiner {
	private HdfsFilesReader reader = null;
	private RecalibratorTable tables = null;
	private Covariate[] covariates = null;

	public RecalibratorTableCombiner(RecalibratorOptions option, SAMFileHeader header) {
		covariates = CovariateUtil.initializeCovariates(option, header);
		tables = new RecalibratorTable(covariates, header.getReadGroups().size());
	}
	
	private class RecalibratorPathFilter implements PathFilter {
		@Override
		public boolean accept(Path path) {
			if (path.getName().startsWith(RecalibratorContextWriter.RECALIBRATOR_TABLE_TAG))
				return true;
			return false;
		}
	}
	
	public static class NonRecalibratorPathFilter implements PathFilter {
		@Override
		public boolean accept(Path path) {
			if (path.getName().startsWith(RecalibratorContextWriter.RECALIBRATOR_TABLE_TAG))
				return false;
			if(path.getName().startsWith("_"))
				return false;
			return true;
		}
	}

	public void combineTable(String path) {
		reader = new HdfsFilesReader();
		reader.traversal(path,new RecalibratorPathFilter());

		String line = null;
		while (reader.hasNext()) {
			line = reader.next();
			tableLineParser(line.split("\t"));
		}
		
		reader.delete();
		reader.clear();
	}

	private void tableLineParser(String[] array) {
		int[] keys = null;
		int length = 2;
		int index = Integer.parseInt(array[0]);
		if (index < 2)
			length += index;
		else
			length = 4;
		keys = new int[length];

		int i;
		for (i = 0; i < length; i++)
			keys[i] = Integer.parseInt(array[i + 1]);

		RecalibratorDatum datum = new RecalibratorDatum(Long.parseLong(array[length + 1]),
				Long.parseLong(array[length + 2]), Double.parseDouble(array[length + 3]));
		
		updateTable(index,keys,datum);
	}
	
	public Covariate[] getCovariates(){
		return this.covariates;
	}

	private void updateTable(int index, int[] keys, RecalibratorDatum datum) {
		final NestedObjectArray<RecalibratorDatum> table = tables.getTable(index);
		final RecalibratorDatum currDatum = table.get(keys);
		if (currDatum == null)
			table.put(datum, keys);
		else {
			if (index == 0)
				currDatum.combine(datum);
			else
				currDatum.increment(datum);
		}
	}
	
	public RecalibratorTable getRecalibratorTable(){
		return this.tables;
	}
}
