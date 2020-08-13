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
package org.bgi.flexlab.gaea.tools.recalibrator;

import htsjdk.samtools.SAMFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.bam.filter.BaseRecalibrationFilter;
import org.bgi.flexlab.gaea.data.structure.reference.BaseAndSNPInformation;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.realigner.RealignerWriter;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil.Consistent;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil.SolidRecallMode;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.CovariateUtil;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTable;
import org.bgi.flexlab.gaea.util.AlignmentUtil;
import org.bgi.flexlab.gaea.util.BaseUtils;
import org.bgi.flexlab.gaea.util.EventType;
import org.bgi.flexlab.gaea.util.NestedObjectArray;

import java.util.ArrayList;

public class RecalibratorEngine {
	private SAMFileHeader mHeader = null;
	private RecalibratorOptions option = new RecalibratorOptions();
	private ReferenceShare chrInfo = null;
	private BaseAndSNPInformation information = null;
	private BaseRecalibrationFilter filter = null;
	private Covariate[] covariates = null;
	private RecalibratorTable recalibratorTables = null;
	private RealignerWriter writer = null;
	private boolean isRealigment = false;

	public RecalibratorEngine(RecalibratorOptions option, ReferenceShare chrInfo, SAMFileHeader mHeader,boolean isRealigment,RealignerWriter writer) {
		this.option = option;
		this.chrInfo = chrInfo;
		this.mHeader = mHeader;
		information = new BaseAndSNPInformation();
		filter = new BaseRecalibrationFilter();
		this.covariates = CovariateUtil.initializeCovariates(option, mHeader);
		recalibratorTables = new RecalibratorTable(this.covariates, mHeader.getReadGroups().size());
		this.isRealigment = isRealigment;
		this.writer = writer;
	}

	private void setWindows(String chrName, int winNum) {
		int start = winNum > 0 ? (winNum - 1) * option.getWindowsSize()+1 : 1;
		int end = (winNum + 2) * option.getWindowsSize();
		information.set(chrInfo, chrName, start, end);
	}

	public void mapReads(ArrayList<GaeaSamRecord> records, Iterable<SamRecordWritable> iterator, String chrName,
			int winNum) {
		setWindows(chrName, winNum);

		if (records != null) {
			for (GaeaSamRecord sam : records) {
				baseQualityStatistics(sam);
			}
		}

		if (iterator == null)
			return;

		for (SamRecordWritable writable : iterator) {
			int readWinNum = writable.get().getAlignmentStart() / option.getWindowsSize();
			GaeaSamRecord sam = new GaeaSamRecord(mHeader, writable.get(), readWinNum == winNum);
			baseQualityStatistics(sam);
			
			if(isRealigment)
				writer.writeRead(sam);
		}
	}

	private Consistent isColorSpaceConsistent(GaeaSamRecord read) {
		return RecalibratorUtil.isColorSpaceConsistent(option.SOLID_NOCALL_STRATEGY, read);
	}

	private void mapRead(GaeaSamRecord read) {
		if (RecalibratorUtil.isSOLiDRead(read) || option.SOLID_RECAL_MODE == SolidRecallMode.DO_NOTHING)
			return;

		int start = read.getAlignmentStart();
		int end = read.getAlignmentEnd();
		int[] readOffsets = AlignmentUtil.readOffsets(read.getCigar(), start, end);

		byte[] quals = read.getBaseQualities();
		byte[] bases = read.getReadBases();

		Consistent consistent = null;
		ReadCovariates rcovariate = null;
		for (int i = start; i <= end; i++) {
			int qpos = readOffsets[i - start];
			if (!information.getSNP(i)) {
				if (consistent == null) {
					consistent = isColorSpaceConsistent(read);
					if (!consistent.isColorSpaceConsistent()) {
						return;
					}
					rcovariate = RecalibratorUtil.computeCovariates(read, covariates);
				}
				if (qpos < 0 || quals[qpos] < option.PRESERVE_QSCORES_LESS_THAN || rcovariate == null
						|| !consistent.isColorSpaceConsistent(qpos, read.getReadNegativeStrandFlag())) {
					continue;
				}
				dataUpdate(qpos, (byte) bases[qpos], quals[qpos], (byte) information.getBase(i), rcovariate);
			}
		}
	}

	public boolean dataUpdate(final int offset, final byte base, final byte quality, final byte refBase,
			final ReadCovariates readCovariates) {
		final boolean isError = !BaseUtils.basesAreEqual(base, refBase);
		final EventType eventType = EventType.SNP;
		final int[] keys = readCovariates.getKeySet(offset, eventType);
		final int eventIndex = eventType.index;

		final NestedObjectArray<RecalibratorDatum> rgRecalTable = recalibratorTables
				.getTable(RecalibratorTable.Type.READ_GROUP_TABLE);

		final RecalibratorDatum rgPreviousDatum = rgRecalTable.get(keys[0], eventIndex);
		final RecalibratorDatum rgThisDatum = RecalibratorDatum.build(quality, isError);
		if (rgPreviousDatum == null)
			rgRecalTable.put(rgThisDatum, keys[0], eventIndex);
		else
			rgPreviousDatum.combine(rgThisDatum);

		final NestedObjectArray<RecalibratorDatum> qualRecalTable = recalibratorTables
				.getTable(RecalibratorTable.Type.QUALITY_SCORE_TABLE);
		final RecalibratorDatum qualPreviousDatum = qualRecalTable.get(keys[0], keys[1], eventIndex);
		if (qualPreviousDatum == null)
			qualRecalTable.put(RecalibratorDatum.build(quality, isError), keys[0], keys[1], eventIndex);
		else
			qualPreviousDatum.increment(isError);

		for (int i = 2; i < covariates.length; i++) {
			if (keys[i] < 0)
				continue;
			final NestedObjectArray<RecalibratorDatum> covRecalTable = recalibratorTables.getTable(i);
			final RecalibratorDatum covPreviousDatum = covRecalTable.get(keys[0], keys[1], keys[i], eventIndex);
			if (covPreviousDatum == null)
				covRecalTable.put(RecalibratorDatum.build(quality, isError), keys[0], keys[1], keys[i], eventIndex);
			else
				covPreviousDatum.increment(isError);
		}
		
		return isError;
	}

	private void baseQualityStatistics(GaeaSamRecord read) {
		if (read.needToOutput()) {
			if (!filter.filter(read, null)) {
				mapRead(read);
			}
		}
	}

	public RecalibratorTable getTables() {
		return this.recalibratorTables;
	}
}
