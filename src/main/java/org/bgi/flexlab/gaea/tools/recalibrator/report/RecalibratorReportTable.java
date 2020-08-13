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

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorDatum;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorUtil;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.quality.QualityQuantizer;
import org.bgi.flexlab.gaea.tools.recalibrator.table.RecalibratorTable;
import org.bgi.flexlab.gaea.util.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RecalibratorReportTable {
	public static final String INVALID_TABLE_NAME_REGEX = "[^a-zA-Z0-9_\\-\\.]";
	public static final int INITITAL_ARRAY_SIZE = 10000;

	private static final String TABLE_HEADER_PREFIX = "#:Table";
	private static final String SEPARATOR = ":";
	private static final String ENDLINE = ":;";
	private static final String ENTER = "\n";
	private static final String TAB = "	";

	private String name;
	private String description;
	private final boolean sortByRowID;

	private Map<Object, Integer> idToIndex;
	private Map<Object, Integer> nameToIndex;
	private List<RecalibratorReportColumn> columns = null;
	private List<Object[]> underlyingData = null;

	public RecalibratorReportTable(final String tableName, final String tableDescription, final int numColumns) {
		this(tableName, tableDescription, numColumns, true);
	}

	public RecalibratorReportTable(final String tableName, final String tableDescription, final int numColumns,
			final boolean sortByRowID) {
		if (!isValidName(tableName)) {
			throw new UserException("Attempted to set a report table name of '" + tableName
					+ "'. Report table names must be purely alphanumeric - no spaces or special characters are allowed.");
		}

		if (!isValidDescription(tableDescription)) {
			throw new UserException("Attempted to set a report table description of '" + tableDescription
					+ "'. Report table descriptions must not contain newlines.");
		}

		this.name = tableName;
		this.description = tableDescription;
		this.sortByRowID = sortByRowID;

		underlyingData = new ArrayList<Object[]>(INITITAL_ARRAY_SIZE);
		columns = new ArrayList<RecalibratorReportColumn>(numColumns);
		nameToIndex = new HashMap<Object, Integer>(numColumns);
		idToIndex = new HashMap<Object, Integer>();
	}

	public RecalibratorReportTable(GaeaFilesReader reader) {
		String[] tableData = null, tableNameData = null;
		if(reader.hasNext()){
			tableData = reader.next().split(SEPARATOR);
		}
		if (reader.hasNext()) {
			tableNameData = reader.next().split(SEPARATOR);
		}

		// parse the header fields
		name = tableNameData[2];
		description = (tableNameData.length <= 3) ? "" : tableNameData[3];
		sortByRowID = false;

		// initialize the data
		final int nColumns = Integer.parseInt(tableData[2]);
		final int nRows = Integer.parseInt(tableData[3]);
		underlyingData = new ArrayList<Object[]>(nRows);
		columns = new ArrayList<RecalibratorReportColumn>(nColumns);
		nameToIndex = new HashMap<Object, Integer>(nColumns);

		// when reading from a file, the row ID mapping is just the index
		idToIndex = new HashMap<Object, Integer>();
		for (int i = 0; i < nRows; i++)
			idToIndex.put(i, i);

		// read the column names
		String columnLine = null;
		if (reader.hasNext())
			columnLine = reader.next();

		final List<Integer> columnStarts = StringUtils.getWordStarts(columnLine);
		final String[] columnNames = StringUtils.splitFixedWidth(columnLine);

		// Put in columns using the format string from the header
		for (int i = 0; i < nColumns; i++) {
			final String format = tableData[4 + i];
			addColumn(columnNames[i], format);
		}

		for (int i = 0; i < nRows; i++) {
			// read a data line
			if (reader.hasNext()) {
				final String dataLine = reader.next();
				final List<String> lineSplits = Arrays.asList(StringUtils.splitFixedWidth(dataLine, columnStarts));

				underlyingData.add(new Object[nColumns]);
				for (int columnIndex = 0; columnIndex < nColumns; columnIndex++) {
					final StandardDataType type = columns.get(columnIndex).getDataType();
					final String columnName = columnNames[columnIndex];
					set(i, columnName, type.Parse(lineSplits.get(columnIndex)));
				}
			}
		}

		if (reader.hasNext())
			reader.next();
	}
	
	public String getTableName(){
		return name;
	}

	private boolean isValidName(String name) {
		Pattern p = Pattern.compile(INVALID_TABLE_NAME_REGEX);
		Matcher m = p.matcher(name);

		return !m.find();
	}

	private boolean isValidDescription(String description) {
		Pattern p = Pattern.compile("\\r|\\n");
		Matcher m = p.matcher(description);

		return !m.find();
	}

	private void expandDataList(final int rowIndex, final boolean updateRowIdMap) {
		int currentSize = underlyingData.size();
		if (rowIndex >= currentSize) {
			final int numNewRows = rowIndex - currentSize + 1;
			for (int i = 0; i < numNewRows; i++) {
				if (updateRowIdMap)
					idToIndex.put(currentSize, currentSize);
				underlyingData.add(new Object[getColumnNumber()]);
				currentSize++;
			}
		}
	}

	public int getColumnNumber() {
		return columns.size();
	}
	
	public int getRowNumber(){
		return underlyingData.size();
	}

	public void addRowID(final String ID) {
		addRowID(ID, false);
	}

	/**
	 * Add a mapping from ID to the index of a new row added to the table.
	 */
	public void addRowID(final String ID, final boolean populateFirstColumn) {
		addRowIDMapping(ID, underlyingData.size(), populateFirstColumn);
	}

	/**
	 * Add a mapping from ID to row index.
	 */
	public void addRowIDMapping(final String ID, final int index) {
		addRowIDMapping(ID, index, false);
	}

	/**
	 * Add a mapping from ID to row index.
	 */
	public void addRowIDMapping(final Object ID, final int index, final boolean populateFirstColumn) {
		expandDataList(index, false);
		idToIndex.put(ID, index);

		if (populateFirstColumn)
			set(index, 0, ID);
	}

	public void set(final Object rowID, final String columnName, final Object value) {
		if (!idToIndex.containsKey(rowID)) {
			idToIndex.put(rowID, underlyingData.size());
			expandDataList(underlyingData.size(), false);
		}
		set(idToIndex.get(rowID), nameToIndex.get(columnName), value);
	}

	public void set(final int rowIndex, final int colIndex, Object value) {
		expandDataList(rowIndex, true);
		verifyEntry(rowIndex, colIndex);
		RecalibratorReportColumn column = columns.get(colIndex);

		if (value == null)
			value = "null";
		else
			value = fixType(value, column);

		if (column.getDataType().equals(StandardDataType.fromObject(value))
				|| column.getDataType().equals(StandardDataType.Unknown)) {
			underlyingData.get(rowIndex)[colIndex] = value;
			column.updateFormatting(value);
		} else {
			throw new UserException(String.format("Tried to add an object of type: %s to a column of type: %s",
					StandardDataType.fromObject(value).name(), column.getDataType().name()));
		}
	}

	private void verifyEntry(final int rowIndex, final int colIndex) {
		if (rowIndex < 0 || colIndex < 0 || colIndex >= getColumnNumber())
			throw new UserException("attempted to access a cell that does not exist in table '" + name + "'");
	}

	private Object fixType(final Object value, final RecalibratorReportColumn column) {
		Object newValue = null;
		if (value instanceof String && !column.getDataType().equals(StandardDataType.String)) {
			newValue = column.getDataType().getNonStringValue(value);
		}

		return (newValue != null) ? newValue : value;
	}

	public void addColumn(String columnName) {
		addColumn(columnName, "");
	}

	public void addColumn(String columnName, String format) {
		nameToIndex.put(columnName, getColumnNumber());
		columns.add(new RecalibratorReportColumn(columnName, format));
	}
	
	public Object get(final int rowIndex, final String columnName) {
        return get(rowIndex, nameToIndex.get(columnName));
    }
	
	public Object get(int rowIndex, int columnIndex) {
        verifyEntry(rowIndex, columnIndex);
        return underlyingData.get(rowIndex)[columnIndex];
    }

	public String toString() {
		StringBuilder tableString = new StringBuilder();

		// print header format
		tableString.append(TABLE_HEADER_PREFIX);
		tableString.append(SEPARATOR);
		tableString.append(String.valueOf(getColumnNumber()));
		tableString.append(SEPARATOR);
		tableString.append(String.valueOf(getRowNumber()));

		for (RecalibratorReportColumn column : columns) {
			tableString.append(SEPARATOR);
			tableString.append(column.getFormat());
		}

		tableString.append(ENDLINE);
		tableString.append(ENTER);

		// print table description
		tableString.append(TABLE_HEADER_PREFIX);
		tableString.append(SEPARATOR);
		tableString.append(name);
		tableString.append(SEPARATOR);
		tableString.append(description);
		tableString.append(ENTER);

		// print column header
		boolean firstColumn = true;
		for (final RecalibratorReportColumn column : columns) {
			if (!firstColumn)
				tableString.append(TAB);
			firstColumn = false;
			tableString.append(String.format(column.getColumnNameFormat(), column.getColumnName()));
		}
		tableString.append(ENTER);

		if (sortByRowID) {
			if (idToIndex.size() != underlyingData.size())
				throw new UserException("There isn't a 1-to-1 mapping from row ID to index!");

			final TreeMap<Object, Integer> sortedMap;
			try {
				sortedMap = new TreeMap<Object, Integer>(idToIndex);
			} catch (ClassCastException e) {
				throw new UserException(
						"Unable to sort the rows based on the row IDs because the ID Objects are of different types");
			}
			for (final Map.Entry<Object, Integer> rowKey : sortedMap.entrySet())
				tableString.append(rowToString(underlyingData.get(rowKey.getValue())));
		} else {
			for (final Object[] row : underlyingData)
				tableString.append(rowToString(row));
		}

		tableString.append(ENTER);

		return tableString.toString();
	}

	private String rowToString(Object[] row) {
		StringBuilder str = new StringBuilder();

		boolean firstColumn = true;
		for (int i = 0; i < row.length; i++) {
			if (!firstColumn)
				str.append(TAB);
			firstColumn = false;

			final Object obj = row[i];
			final String value;

			final RecalibratorReportColumn info = columns.get(i);

			if (obj == null)
				value = "null";
			else if (info.getDataType().equals(StandardDataType.Unknown)
					&& (obj instanceof Double || obj instanceof Float))
				value = String.format("%.8f", obj);
			else
				value = String.format(info.getFormat(), obj);

			str.append(String.format(info.getColumnValueFormat(), value));
		}

		str.append(ENTER);

		return str.toString();
	}

	public static RecalibratorReportTable reportTableBuilder(RecalibratorOptions option, String covariateNames) {
		RecalibratorReportTable argsTable = new RecalibratorReportTable(RecalibratorUtil.ARGUMENT_TABLE_NAME,
				"Recalibration argument collection values used in this run", 2);
		argsTable.addColumn("Argument");
		argsTable.addColumn(RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME);
		argsTable.addRowID("covariate", true);
		argsTable.set("covariate", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, covariateNames);
		argsTable.addRowID("no_standard_covs", true);
		argsTable.set("no_standard_covs", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.DO_NOT_USE_STANDARD_COVARIATES);
		argsTable.addRowID("run_without_dbsnp", true);
		argsTable.addRowID("solid_recal_mode", true);
		argsTable.set("solid_recal_mode", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.SOLID_RECAL_MODE);
		argsTable.addRowID("solid_nocall_strategy", true);
		argsTable.set("solid_nocall_strategy", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.SOLID_NOCALL_STRATEGY);
		argsTable.addRowID("mismatches_context_size", true);
		argsTable.set("mismatches_context_size", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.MISMATCHES_CONTEXT_SIZE);
		argsTable.addRowID("indels_context_size", true);
		argsTable.set("indels_context_size", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.INDELS_CONTEXT_SIZE);
		argsTable.addRowID("mismatches_default_quality", true);
		argsTable.set("mismatches_default_quality", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.MISMATCHES_DEFAULT_QUALITY);
		argsTable.addRowID("insertions_default_quality", true);
		argsTable.set("insertions_default_quality", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.INSERTIONS_DEFAULT_QUALITY);
		argsTable.addRowID("low_quality_tail", true);
		argsTable.set("low_quality_tail", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.LOW_QUALITY_TAIL);
		argsTable.addRowID("default_platform", true);
		argsTable.set("default_platform", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.DEFAULT_PLATFORM);
		argsTable.addRowID("force_platform", true);
		argsTable.set("force_platform", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.FORCE_PLATFORM);
		argsTable.addRowID("quantizing_levels", true);
		argsTable.set("quantizing_levels", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.QUANTIZING_LEVELS);
		argsTable.addRowID("keep_intermediate_files", true);
		argsTable.set("keep_intermediate_files", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.KEEP_INTERMEDIATE_FILES);
		argsTable.addRowID("no_plots", true);
		argsTable.set("no_plots", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME, option.NO_PLOTS);
		argsTable.addRowID("recalibration_report", true);
		argsTable.addRowID("binary_tag_name", true);
		argsTable.set("binary_tag_name", RecalibratorUtil.ARGUMENT_VALUE_COLUMN_NAME,
				option.BINARY_TAG_NAME == null ? "null" : option.BINARY_TAG_NAME);
		return argsTable;
	}

	public static RecalibratorReportTable reportTableBuilder(RecalibratorTable tables, int nLevels) {
		final Long[] qualityHistogram = new Long[QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE + 1];

		for (int i = 0; i < qualityHistogram.length; i++)
			qualityHistogram[i] = 0L;

		final NestedObjectArray<RecalibratorDatum> qualTable = tables
				.getTable(RecalibratorTable.Type.QUALITY_SCORE_TABLE);

		for (final RecalibratorDatum value : qualTable.getAllValues()) {
			final RecalibratorDatum datum = value;
			final int empiricalQual = MathUtils.fastRound(datum.getEmpiricalQuality());

			qualityHistogram[empiricalQual] += datum.getBasesNumber();

		}
		List<Long> empiricalCounts = Arrays.asList(qualityHistogram);
		QualityQuantizer quantizer = new QualityQuantizer(empiricalCounts, nLevels,
				QualityUtils.MINIMUM_USABLE_QUALITY_SCORE);
		List<Byte> quantizedQualities = quantizer.getIntervals();

		return reportTableBuilder(empiricalCounts, quantizedQualities);
	}

	private static RecalibratorReportTable reportTableBuilder(List<Long> empiricalCounts, List<Byte> quantizedCounts) {
		RecalibratorReportTable quantizedTable = new RecalibratorReportTable(RecalibratorUtil.QUANTIZED_TABLE_NAME,
				"Quality quantization map", 3);
		quantizedTable.addColumn(RecalibratorUtil.QUANTIZED_SCORE_COLUMN_NAME);
		quantizedTable.addColumn(RecalibratorUtil.QUANTIZED_COUNT_COLUMN_NAME);
		quantizedTable.addColumn(RecalibratorUtil.QUANTIZED_VALUE_COLUMN_NAME);
		for (int qual = 0; qual <= QualityUtils.MAXIMUM_USABLE_QUALITY_SCORE; qual++) {
			quantizedTable.set(qual, RecalibratorUtil.QUANTIZED_SCORE_COLUMN_NAME, qual);
			quantizedTable.set(qual, RecalibratorUtil.QUANTIZED_COUNT_COLUMN_NAME, empiricalCounts.get(qual));
			quantizedTable.set(qual, RecalibratorUtil.QUANTIZED_VALUE_COLUMN_NAME, quantizedCounts.get(qual));
		}
		return quantizedTable;
	}

	private static String parseCovariateName(final Covariate covariate) {
		return covariate.getClass().getSimpleName().split("Covariate")[0];
	}

	private static final Pair<String, String> covariateValue = new Pair<String, String>(
			RecalibratorUtil.COVARIATE_VALUE_COLUMN_NAME, "%s");
	private static final Pair<String, String> covariateName = new Pair<String, String>(
			RecalibratorUtil.COVARIATE_NAME_COLUMN_NAME, "%s");
	private static final Pair<String, String> eventType = new Pair<String, String>(
			RecalibratorUtil.EVENT_TYPE_COLUMN_NAME, "%s");
	private static final Pair<String, String> empiricalQuality = new Pair<String, String>(
			RecalibratorUtil.EMPIRICAL_QUALITY_COLUMN_NAME, "%.4f");
	private static final Pair<String, String> estimatedQReported = new Pair<String, String>(
			RecalibratorUtil.ESTIMATED_Q_REPORTED_COLUMN_NAME, "%.4f");
	private static final Pair<String, String> nObservations = new Pair<String, String>(
			RecalibratorUtil.NUMBER_OBSERVATIONS_COLUMN_NAME, "%d");
	private static final Pair<String, String> nErrors = new Pair<String, String>(
			RecalibratorUtil.NUMBER_ERRORS_COLUMN_NAME, "%d");

	public static List<RecalibratorReportTable> reportTableBuilder(RecalibratorTable tables, Covariate[] covariates) {
		List<RecalibratorReportTable> tableList = new ArrayList<RecalibratorReportTable>();

		int reportTableIndex = 0;
		int rowIndex = 0;

		final Map<Covariate, String> covariateNameMap = new HashMap<Covariate, String>(covariates.length);
		for (final Covariate covariate : covariates)
			covariateNameMap.put(covariate, parseCovariateName(covariate));

		for (int tableIndex = 0; tableIndex < tables.length(); tableIndex++) {
			final ArrayList<Pair<String, String>> columnNames = new ArrayList<Pair<String, String>>();
			columnNames.add(new Pair<String, String>(covariateNameMap.get(covariates[0]), "%s"));
			if (tableIndex != RecalibratorTable.Type.READ_GROUP_TABLE.index) {
				columnNames.add(new Pair<String, String>(covariateNameMap.get(covariates[1]), "%s"));
				if (tableIndex >= RecalibratorTable.Type.OPTIONAL_COVARIATE_TABLES_START.index) {
					columnNames.add(covariateValue);
					columnNames.add(covariateName);
				}
			}

			columnNames.add(eventType);
			columnNames.add(empiricalQuality);
			if (tableIndex == RecalibratorTable.Type.READ_GROUP_TABLE.index)
				columnNames.add(estimatedQReported);
			columnNames.add(nObservations);
			columnNames.add(nErrors);

			final RecalibratorReportTable reportTable;
			if (tableIndex <= RecalibratorTable.Type.OPTIONAL_COVARIATE_TABLES_START.index) {
				reportTable = new RecalibratorReportTable(RecalibratorUtil.RECALIBRATOR_TABLE_NAME[reportTableIndex++], "", columnNames.size());
				for (final Pair<String, String> columnName : columnNames)
					reportTable.addColumn(columnName.getFirst(), columnName.getSecond());
				rowIndex = 0;
			} else {
				reportTable = tableList.get(RecalibratorTable.Type.OPTIONAL_COVARIATE_TABLES_START.index);
			}

			final NestedObjectArray<RecalibratorDatum> table = tables.getTable(tableIndex);
			for (final NestedObjectArray.Leave row : table.getAllLeaves()) {
				final RecalibratorDatum datum = (RecalibratorDatum) row.value;
				final int[] keys = row.keys;

				int columnIndex = 0;
				int keyIndex = 0;
				reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(),
						covariates[0].formatKey(keys[keyIndex++]));
				if (tableIndex != RecalibratorTable.Type.READ_GROUP_TABLE.index) {
					reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(),
							covariates[1].formatKey(keys[keyIndex++]));
					if (tableIndex >= RecalibratorTable.Type.OPTIONAL_COVARIATE_TABLES_START.index) {
						final Covariate covariate = covariates[tableIndex];

						reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(),
								covariate.formatKey(keys[keyIndex++]));
						reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(),
								covariateNameMap.get(covariate));
					}
				}

				final EventType event = EventType.eventFrom(keys[keyIndex]);
				reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(), event.toString());

				reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(), datum.getEmpiricalQuality());
				if (tableIndex == RecalibratorTable.Type.READ_GROUP_TABLE.index)
					reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(), datum.getEstimatedQuality());
				reportTable.set(rowIndex, columnNames.get(columnIndex++).getFirst(), datum.getBasesNumber());
				reportTable.set(rowIndex, columnNames.get(columnIndex).getFirst(), datum.getMismatchNumber());

				rowIndex++;
			}
			tableList.add(reportTable);
		}

		return tableList;
	}
}
