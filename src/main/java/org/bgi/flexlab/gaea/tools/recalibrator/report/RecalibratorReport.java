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
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.util.HdfsFilesReader;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamTag;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.ReadCovariates;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorDatum;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.CovariateUtil;
import org.bgi.flexlab.gaea.tools.recalibrator.quality.QualityQuantizer;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTable;
import org.bgi.flexlab.gaea.util.*;

import java.util.*;

public class RecalibratorReport {
	public static final int MAXIMUM_RECALIBRATED_READ_LENGTH = 2000;

	private HashMap<String, RecalibratorReportTable> tables = null;
	private RecalibratorOptions option = new RecalibratorOptions();
	private HashMap<String, Integer> optionalIndex = new HashMap<String, Integer>();
	private RecalibratorTable recalTable = null;
	private Covariate[] covariates = null;

	private List<Byte> qualities = null;
	private List<Long> qualityCounts = null;
	private int quanLevels;

	private final boolean diableIndelQuality;
	private final boolean emitOriginQuality;
	private final int preserveQualityLessThan;

	private final ReadCovariates readCovariates;

	public RecalibratorReport(String input, SAMFileHeader header, int quanLevels, int preserveQualityLessThan,
			boolean disableIndelQuality, boolean emitOriginQuality) {
		initialize(input, header);

		qualityInitialized(quanLevels);

		readCovariates = new ReadCovariates(MAXIMUM_RECALIBRATED_READ_LENGTH, covariates.length);

		this.diableIndelQuality = disableIndelQuality;
		this.preserveQualityLessThan = preserveQualityLessThan;
		this.emitOriginQuality = emitOriginQuality;
	}

	public RecalibratorReport(String input, SAMFileHeader header, int quanLevels, int preserveQualityLessThan) {
		this(input, header, quanLevels, preserveQualityLessThan, false, false);
	}

	private void initialize(String input, SAMFileHeader header) {
		tables = new HashMap<String, RecalibratorReportTable>();
		addTables(input);
		option.parse(tables.get(RecalibratorUtil.ARGUMENT_TABLE_NAME));

		covariates = CovariateUtil.initializeCovariates(option, header);
		for (int i = 2; i < covariates.length; i++) {
			String covName = covariates[i].getClass().getSimpleName().split("Covariate")[0];
			optionalIndex.put(covName, i - 2);
		}

		recalTable = new RecalibratorTable(covariates,
				readGroupSize(tables.get(RecalibratorUtil.RECALIBRATOR_TABLE_NAME[0])));
		// read group table parse
		readGroupParser(tables.get(RecalibratorUtil.RECALIBRATOR_TABLE_NAME[0]),
				recalTable.getTable(RecalibratorTable.Type.READ_GROUP_TABLE));
		
		// quality group table parse
		qualityScoreParser(tables.get(RecalibratorUtil.RECALIBRATOR_TABLE_NAME[1]),
				recalTable.getTable(RecalibratorTable.Type.QUALITY_SCORE_TABLE));
		// covariate tables parse
		covariateParser(tables.get(RecalibratorUtil.RECALIBRATOR_TABLE_NAME[2]), recalTable);
	}

	private void qualityInitialized(int levels) {
		this.quanLevels = levels;

		RecalibratorReportTable qualTable = tables.get(RecalibratorUtil.QUANTIZED_TABLE_NAME);

		final Byte[] quals = new Byte[QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE + 1];
		final Long[] counts = new Long[QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE + 1];

		for (int i = 0; i < qualTable.getRowNumber(); i++) {
			final byte originalQual = (byte) i;
			final Object quantizedObject = qualTable.get(i, RecalibratorUtil.QUANTIZED_VALUE_COLUMN_NAME);
			final Object countObject = qualTable.get(i, RecalibratorUtil.QUANTIZED_COUNT_COLUMN_NAME);
			final byte quantizedQual = Byte.parseByte(quantizedObject.toString());
			final long quantizedCount = Long.parseLong(countObject.toString());
			quals[originalQual] = quantizedQual;
			counts[originalQual] = quantizedCount;
		}

		qualities = Arrays.asList(quals);
		qualityCounts = Arrays.asList(counts);

		int tmpLevel = calculateLevels(qualities);

		if (levels == 0) {
			this.quanLevels = QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE;
			for (int i = 0; i < this.quanLevels; i++)
				qualities.set(i, (byte) i);
		} else if (levels > 0 && levels != tmpLevel) {
			QualityQuantizer quantizer = new QualityQuantizer(qualityCounts, levels,
					QualityUtils.MINIMUM_USABLE_QUALITY_SCORE);
			qualities = quantizer.getIntervals();
		}
	}

	private int calculateLevels(List<Byte> quals) {
		byte lastByte = -1;
		int levels = 0;
		for (byte q : quals) {
			if (q != lastByte) {
				levels++;
				lastByte = q;
			}
		}
		return levels;
	}

	public int readGroupSize(RecalibratorReportTable table) {
		Set<String> readGroup = new HashSet<String>();

		for (int i = 0; i < table.getRowNumber(); i++)
			readGroup.add(table.get(i, RecalibratorUtil.READGROUP_COLUMN_NAME).toString());
		return readGroup.size();
	}

	public void addTables(String input) {
		HdfsFilesReader reader = new HdfsFilesReader();
		reader.traversal(input);

		String header = null;
		if (reader.hasNext()) {
			header = reader.next();
		}
		if (header == null)
			throw new RuntimeException("heade is null");

		int ntables = Integer.parseInt(header.split(":")[2]);

		for (int i = 0; i < ntables; i++) {
			RecalibratorReportTable table = new RecalibratorReportTable(reader);
			tables.put(table.getTableName(), table);
		}
	}

	private void readGroupParser(RecalibratorReportTable table, NestedObjectArray<RecalibratorDatum> nestedArray) {
		final int[] rgArray = new int[2];

		for (int i = 0; i < table.getRowNumber(); i++) {
			final Object rg = table.get(i, RecalibratorUtil.READGROUP_COLUMN_NAME);
			rgArray[0] = covariates[0].keyFromValue(rg);
			final EventType event = EventType.eventFrom((String) table.get(i, RecalibratorUtil.EVENT_TYPE_COLUMN_NAME));
			rgArray[1] = event.index;

			nestedArray.put(RecalibratorDatum.build(table, i, true), rgArray);
		}
	}

	private void qualityScoreParser(RecalibratorReportTable table, NestedObjectArray<RecalibratorDatum> nestedArray) {
		final int[] qualArray = new int[3];

		for (int i = 0; i < table.getRowNumber(); i++) {
			final Object rg = table.get(i, RecalibratorUtil.READGROUP_COLUMN_NAME);
			qualArray[0] = covariates[0].keyFromValue(rg);
			final Object qual = table.get(i, RecalibratorUtil.QUANTIZED_SCORE_COLUMN_NAME);
			qualArray[1] = covariates[1].keyFromValue(qual);
			final EventType event = EventType.eventFrom((String) table.get(i, RecalibratorUtil.EVENT_TYPE_COLUMN_NAME));
			qualArray[2] = event.index;

			nestedArray.put(RecalibratorDatum.build(table, i, false), qualArray);
		}
	}

	private void covariateParser(RecalibratorReportTable table, RecalibratorTable recalTable) {
		final int[] covArray = new int[4];

		for (int i = 0; i < table.getRowNumber(); i++) {
			final Object rg = table.get(i, RecalibratorUtil.READGROUP_COLUMN_NAME);
			covArray[0] = covariates[0].keyFromValue(rg);
			final Object qual = table.get(i, RecalibratorUtil.QUANTIZED_SCORE_COLUMN_NAME);
			covArray[1] = covariates[1].keyFromValue(qual);

			final String covName = (String) table.get(i, RecalibratorUtil.COVARIATE_NAME_COLUMN_NAME);
			final int covIndex = optionalIndex.get(covName);
			final Object covValue = table.get(i, RecalibratorUtil.COVARIATE_VALUE_COLUMN_NAME);
			covArray[2] = covariates[RecalibratorTable.Type.OPTIONAL_COVARIATE_TABLES_START.index + covIndex]
					.keyFromValue(covValue);

			final EventType event = EventType.eventFrom((String) table.get(i, RecalibratorUtil.EVENT_TYPE_COLUMN_NAME));
			covArray[3] = event.index;

			recalTable.getTable(RecalibratorTable.Type.OPTIONAL_COVARIATE_TABLES_START.index + covIndex)
					.put(RecalibratorDatum.build(table, i, false), covArray);
		}
	}

	public void readRecalibrator(GaeaSamRecord read) {
		if (this.emitOriginQuality && read.getAttribute(GaeaSamTag.OQ.name()) == null) {
			try {
				read.setAttribute(GaeaSamTag.OQ.name(), SAMUtils.phredToFastq(read.getBaseQualities()));
			} catch (IllegalArgumentException e) {
				throw new UserException.MalformedBAM(read, "illegal base quality encountered; " + e.getMessage());
			}
		}

		RecalibratorUtil.computeCovariates(read, covariates, readCovariates);
		for (final EventType errorModel : EventType.values()) {
			if (errorModel != EventType.SNP) {
				continue; // current only support SNP event
			}

			if (this.diableIndelQuality && errorModel != EventType.SNP) {
				read.setBaseQualities(null, errorModel);
				continue;
			}

			final byte[] qualities = read.getBaseQualities(errorModel);
			final int[][] fullReadKeySet = readCovariates.getKeySet(errorModel);

			final int readLength = read.getReadLength();
			for (int offset = 0; offset < readLength; offset++) {

				final byte originalQualityScore = qualities[offset];

				if (originalQualityScore >= this.preserveQualityLessThan) {
					final int[] keySet = fullReadKeySet[offset];
					final byte recalibratedQualityScore = performSequentialQualityCalculation(keySet, errorModel);
					qualities[offset] = recalibratedQualityScore;
				}
			}
			
			read.setBaseQualities(qualities, errorModel);
		}
	}

	private byte performSequentialQualityCalculation(final int[] key, final EventType errorModel) {
		final byte qualFromRead = (byte) (long) key[1];
		final double globalDeltaQ = calculateGlobalDeltaQ(recalTable.getTable(RecalibratorTable.Type.READ_GROUP_TABLE),
				key, errorModel);
		final double deltaQReported = calculateDeltaQReported(
				recalTable.getTable(RecalibratorTable.Type.QUALITY_SCORE_TABLE), key, errorModel, globalDeltaQ,
				qualFromRead);
		final double deltaQCovariates = calculateDeltaQCovariates(recalTable, key, errorModel, globalDeltaQ,
				deltaQReported, qualFromRead);

		double recalibratedQual = qualFromRead + globalDeltaQ + deltaQReported + deltaQCovariates;
		recalibratedQual = QualityUtils.boundQuality(MathUtils.fastRound(recalibratedQual),
				QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE);

		return this.qualities.get((int) recalibratedQual);
	}

	private double calculateGlobalDeltaQ(final NestedObjectArray<RecalibratorDatum> table, final int[] key,
			final EventType errorModel) {
		double result = 0.0;

		final RecalibratorDatum empiricalQualRG = table.get(key[0], errorModel.index);
		if (empiricalQualRG != null) {
			final double globalDeltaQEmpirical = empiricalQualRG.getEmpiricalQuality();
			final double aggregrateQReported = empiricalQualRG.getEstimatedQuality();
			result = globalDeltaQEmpirical - aggregrateQReported;
		}

		return result;
	}

	private double calculateDeltaQReported(final NestedObjectArray<RecalibratorDatum> table, final int[] key,
			final EventType errorModel, final double globalDeltaQ, final byte qualFromRead) {
		double result = 0.0;

		final RecalibratorDatum empiricalQualQS = table.get(key[0], key[1], errorModel.index);
		if (empiricalQualQS != null) {
			final double deltaQReportedEmpirical = empiricalQualQS.getEmpiricalQuality();
			result = deltaQReportedEmpirical - qualFromRead - globalDeltaQ;
		}

		return result;
	}

	private double calculateDeltaQCovariates(final RecalibratorTable recalibrationTables, final int[] key,
			final EventType errorModel, final double globalDeltaQ, final double deltaQReported,
			final byte qualFromRead) {
		double result = 0.0;

		// for all optional covariates
		for (int i = 2; i < covariates.length; i++) {
			if (key[i] < 0)
				continue;

			final RecalibratorDatum empiricalQualCO = recalibrationTables.getTable(i).get(key[0], key[1], key[i],
					errorModel.index);
			if (empiricalQualCO != null) {
				final double deltaQCovariateEmpirical = empiricalQualCO.getEmpiricalQuality();
				result += (deltaQCovariateEmpirical - qualFromRead - (globalDeltaQ + deltaQReported));
			}
		}
		return result;
	}

	public void clear() {
		if (tables != null)
			tables.clear();
		optionalIndex.clear();
	}
}
