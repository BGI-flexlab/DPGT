package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.report;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util.BucketUtils;

public final class Report {
	public static final String RECAL_FILE = "input covariates table file for base quality score recalibration";
    public static final String REPORT_HEADER_PREFIX = "#:Report.";
    public static final ReportVersion LATEST_REPORT_VERSION = ReportVersion.V1_1;
    private static final String SEPARATOR = ":";
    private ReportVersion version = LATEST_REPORT_VERSION;
    
    public static final String READGROUP_COLUMN_NAME = "ReadGroup";
    public static final String READGROUP_REPORT_TABLE_TITLE = "RecalTable0";

    private final NavigableMap<String, ReportTable> tables = new TreeMap<>();

    /**
     * Create a new, empty Report.
     */
    public Report() {
    }

    /**
     * Create a new Report with the contents of a Report on disk.
     *
     * @param filename the path to the file to load
     */
    public Report(String filename) {
        this(BucketUtils.openFile(filename));
    }

    /**
     * Create a new Report with the contents of a Report on disk.
     *
     * @param file the file to load
     */
    public Report(File file) {
        this(file.getPath());
    }

    public Report(InputStream in){
        loadReport(new InputStreamReader(in));
    }

    /**
     * Create a new report from report tables
     * @param tables Any number of tables that you want to add to the report
     */
    public Report(ReportTable... tables) {
        for( ReportTable table: tables)
            addTable(table);
    }

    /**
     * Gets the unique read groups in the table
     *
     * @return the unique read groups
     */
    public SortedSet<String> getReadGroups() {
        final ReportTable reportTable = getTable(READGROUP_REPORT_TABLE_TITLE);
        final SortedSet<String> readGroups = new TreeSet<>();
        for ( int i = 0; i < reportTable.getNumRows(); i++ ) {
            readGroups.add(reportTable.get(i, READGROUP_COLUMN_NAME).toString());
        }
        return readGroups;
    }

    /**
     * Load a Report from a {@link Reader}
     *
     * @param in the reader to load from
     */
    private void loadReport(Reader in) {
        BufferedReader reader = new BufferedReader(in);
        String reportHeader;
        try {
            reportHeader = reader.readLine();
        } catch (IOException e) {
            throw new UserException("Could not read " + RECAL_FILE, e);
        }

        if ( reportHeader == null ) {
            throw new UserException(RECAL_FILE + " is empty.");
        }

        // Read the first line for the version and number of tables.
        version = ReportVersion.fromHeader(reportHeader);
        if (version.equals(ReportVersion.V0_1) ||
                version.equals(ReportVersion.V0_2))
            throw new UserException("No longer supports reading legacy Reports. Please use v1.0 or newer.");

        int nTables = Integer.parseInt(reportHeader.split(":")[2]);

        // Read each table according ot the number of tables
        for (int i = 0; i < nTables; i++) {
            addTable(new ReportTable(reader, version));
        }
    }


    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param numColumns       the number of columns in this table
     */
    public void addTable(final String tableName, final String tableDescription, final int numColumns) {
        addTable(tableName, tableDescription, numColumns, ReportTable.Sorting.DO_NOT_SORT);
    }

    /**
     * Add a new, empty table to the report
     *
     * @param tableName        the name of the table
     * @param tableDescription the description of the table
     * @param numColumns       the number of columns in this table
     * @param sortingWay       way to sort table
     */
    public void addTable(final String tableName, final String tableDescription, final int numColumns, final ReportTable.Sorting sortingWay) {
    	ReportTable table = new ReportTable(tableName, tableDescription, numColumns, sortingWay);
        tables.put(tableName, table);
    }

    /**
     * Adds a table, empty or populated, to the report
     *
     * @param table the table to add
     */
    public void addTable(ReportTable table) {
        tables.put(table.getTableName(), table);
    }

    public void addTables(List<ReportTable> ReportTableV2s) {
        for ( ReportTable table : ReportTableV2s )
            addTable(table);
    }

    /**
     * Return a table with a given name
     *
     * @param tableName the name of the table
     * @return the table object
     */
    public ReportTable getTable(String tableName) {
    	ReportTable table = tables.get(tableName);
        if (table == null)
            throw new UserException("Table is not in Report: " + tableName);
        return table;
    }

    /**
     * Print all tables contained within this container to a PrintStream
     *
     * @param out the PrintStream to which the tables should be written
     */
    public void print(PrintStream out) {
        out.println(REPORT_HEADER_PREFIX + getVersion() + SEPARATOR + getTables().size());
        for (ReportTable table : tables.values()) {
            table.write(out);
        }
    }

    /**
     * Print all tables contained within this container to a PrintStream
     *
     * @param out the PrintStream to which the tables should be written
     */
    public void print(PrintStream out, ReportTable.Sorting sortingWay) {
        out.println(REPORT_HEADER_PREFIX + getVersion() + SEPARATOR + getTables().size());
        for (ReportTable table : tables.values()) {
            table.write(out, sortingWay);
        }
    }



    public Collection<ReportTable> getTables() {
        return tables.values();
    }

    /**
     * This is the main function is charge of gathering the reports. It checks that the reports are compatible and then
     * calls the table gathering functions.
     *
     * @param input another Report of the same format
     */
    public void concat(Report input) {

        if ( !isSameFormat(input) ) {
            throw new UserException("Failed to combine Report, format doesn't match!");
        }

        for ( Map.Entry<String, ReportTable> table : tables.entrySet() ) {
            table.getValue().concat(input.getTable(table.getKey()));
        }
    }

    public ReportVersion getVersion() {
        return version;
    }

    /**
     * Returns whether or not the two reports have the same format, from columns, to tables, to reports, and everything
     * in between. This does not check if the data inside is the same. This is the check to see if the two reports are
     * gatherable or reduceable.
     *
     * @param report another report
     * @return true if the the reports are gatherable
     */
    public boolean isSameFormat(Report report) {
        if (!version.equals(report.version)) {
            return false;
        }
        if (!tables.keySet().equals(report.tables.keySet())) {
            return false;
        }
        for (String tableName : tables.keySet()) {
            if (!getTable(tableName).isSameFormat(report.getTable(tableName)))
                return false;
        }
        return true;
    }

    /**
     * Checks that the reports are exactly the same.
     *
     * @param report another report
     * @return true if all field in the reports, tables, and columns are equal.
     */
    public boolean equals(Report report) {
        if (!version.equals(report.version)) {
            return false;
        }
        if (!tables.keySet().equals(report.tables.keySet())) {
            return false;
        }
        for (String tableName : tables.keySet()) {
            if (!getTable(tableName).equals(report.getTable(tableName)))
                return false;
        }
        return true;
    }

    /**
     * The constructor for a simplified Report. Simplified report are designed for reports that do not need
     * the advanced functionality of a full Report.
     * <p/>
     * A simple Report consists of:
     * <p/>
     * - A single table
     * - No primary key ( it is hidden )
     * <p/>
     * Optional:
     * - Only untyped columns. As long as the data is an Object, it will be accepted.
     * - Default column values being empty strings.
     * <p/>
     * Limitations:
     * <p/>
     * - A simple report cannot contain multiple tables.
     * - It cannot contain typed columns, which prevents arithmetic gathering.
     *
     * @param tableName The name of your simplereport table
     * @param columns   The names of the columns in your table
     * @return a simplified report
     */
    public static Report newSimpleReport(final String tableName, ReportTable.Sorting sorting, final String... columns) {
        return newSimpleReportWithDescription(tableName, "A simplified table report", sorting, columns);
    }

    /**
     * @see #newSimpleReport(String, ReportTable.Sorting, String...) but with a customized description
     */
    public static Report newSimpleReportWithDescription(final String tableName, final String desc, ReportTable.Sorting sorting, final String... columns) {
    	ReportTable table = new ReportTable(tableName, desc, columns.length, sorting);

        for (String column : columns) {
            table.addColumn(column, "");
        }

        Report output = new Report();
        output.addTable(table);

        return output;
    }

    /**
     * This method provides an efficient way to populate a simplified report. This method will only work on reports
     * that qualify as simplified reports. See the newSimpleReport() constructor for more information.
     *
     * @param values     the row of data to be added to the table.
     *               Note: the number of arguments must match the columns in the table.
     */
    public void addRow(final Object... values) {
        // Must be a simple report
        if ( tables.size() != 1 )
            throw new UserException("Cannot write a row to a complex Report");

        ReportTable table = tables.firstEntry().getValue();
        if ( table.getNumColumns() != values.length )
            throw new UserException("The number of arguments in writeRow (" + values.length + ") must match the number of columns in the table (" + table.getNumColumns() + ")" );

        final int rowIndex = table.getNumRows();
        for ( int i = 0; i < values.length; i++ )
            table.set(rowIndex, i, values[i]);
    }
}
