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
package org.bgi.flexlab.gaea.tools.recalibrator.report;

import htsjdk.samtools.SAMFileHeader;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTable;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTableCombiner;

import java.util.ArrayList;
import java.util.List;

public class RecalibratorReportTableEngine {
	private RecalibratorOptions option = null;
	private SAMFileHeader header = null;
	private List<RecalibratorReportTable> reportTables = null;
	private RecalibratorReportWriter writer = null;

	public RecalibratorReportTableEngine(RecalibratorOptions option, SAMFileHeader header,
			RecalibratorReportWriter writer) {
		this.header = header;
		this.option = option;
		this.writer = writer;
	}

	public void writeReportTable(String input) {
		getReportTables(input);
		print();
		clear();
	}

	private void getReportTables(String input) {
		RecalibratorTableCombiner combiner = new RecalibratorTableCombiner(option, header);
		combiner.combineTable(input);
		RecalibratorTable table = combiner.getRecalibratorTable();

		reportTables = new ArrayList<RecalibratorReportTable>();

		reportTables.add(RecalibratorReportTable.reportTableBuilder(option, covariateNames(combiner.getCovariates())));
		reportTables.add(RecalibratorReportTable.reportTableBuilder(table, option.QUANTIZING_LEVELS));
		reportTables.addAll(RecalibratorReportTable.reportTableBuilder(table, combiner.getCovariates()));
	}

	private void print() {
		writer.write("#:Report.V1.0:" + reportTables.size() + "\n");
		for (RecalibratorReportTable table : reportTables) {
			writer.write(table);
		}
	}

	private void clear() {
		reportTables.clear();
		writer.close();
	}

	private String covariateNames(Covariate[] covariates) {
		boolean first = true;
		StringBuilder sb = new StringBuilder();
		for (final Covariate cov : covariates) {
			if (first) {
				sb.append(cov.getClass().getSimpleName());
				first = false;
			} else {
				sb.append("," + cov.getClass().getSimpleName());
			}
		}

		return sb.toString();
	}
}
