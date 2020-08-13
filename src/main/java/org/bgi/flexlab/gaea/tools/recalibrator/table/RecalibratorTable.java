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
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorDatum;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.CovariateUtil;
import org.bgi.flexlab.gaea.util.EventType;
import org.bgi.flexlab.gaea.util.NestedObjectArray;
import org.bgi.flexlab.gaea.util.NestedObjectArray.Leave;

import java.util.ArrayList;
import java.util.List;

public class RecalibratorTable {
	public enum Type {
		READ_GROUP_TABLE(0), QUALITY_SCORE_TABLE(1), OPTIONAL_COVARIATE_TABLES_START(2);

		public final int index;

		private Type(final int index) {
			this.index = index;
		}
	}

	@SuppressWarnings("rawtypes")
	private NestedObjectArray[] tables = null;

	public RecalibratorTable(final Covariate[] covariates, int readGroupNumber) {
		tables = new NestedObjectArray[covariates.length];

		int maxQualityScore = covariates[Type.QUALITY_SCORE_TABLE.index].maximumKeyValue() + 1;
		int eventSize = EventType.values().length;

		tables[Type.READ_GROUP_TABLE.index] = new NestedObjectArray<RecalibratorDatum>(readGroupNumber, eventSize);
		tables[Type.QUALITY_SCORE_TABLE.index] = new NestedObjectArray<RecalibratorDatum>(readGroupNumber,
				maxQualityScore, eventSize);
		for (int i = Type.OPTIONAL_COVARIATE_TABLES_START.index; i < covariates.length; i++)
			tables[i] = new NestedObjectArray<RecalibratorDatum>(readGroupNumber, maxQualityScore,
					covariates[i].maximumKeyValue() + 1, eventSize);
	}

	public static RecalibratorTable build(RecalibratorOptions option, SAMFileHeader header) {
		RecalibratorTable recalibratorTables = new RecalibratorTable(CovariateUtil.initializeCovariates(option, header),
				header.getReadGroups().size());

		return recalibratorTables;
	}

	@SuppressWarnings("unchecked")
	public NestedObjectArray<RecalibratorDatum> getTable(int index) {
		return (NestedObjectArray<RecalibratorDatum>) tables[index];
	}

	public NestedObjectArray<RecalibratorDatum> getTable(Type type) {
		return getTable(type.index);
	}

	public int length() {
		return tables.length;
	}

	public ArrayList<String> valueStrings() {
		ArrayList<String> arrays = new ArrayList<String>();
		for (int i = 0; i < tables.length; i++) {
			@SuppressWarnings("unchecked")
			NestedObjectArray<RecalibratorDatum> table = tables[i];
			List<Leave> leaves = table.getAllLeaves();
			for (Leave leave : leaves) {
				arrays.add(leave.toString(i));
			}
		}
		return arrays;
	}
}
