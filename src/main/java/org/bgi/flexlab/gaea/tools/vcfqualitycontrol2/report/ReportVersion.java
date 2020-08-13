package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.report;

import org.bgi.flexlab.gaea.data.exception.UserException;

public enum ReportVersion {
    /**
     * Differences between other versions:
     * - Does not allow spaces in cells.
     * - Mostly fixed width but has a bug where the string width of floating point
     * values was not measured correctly leading to columns that aren't aligned
     */
    V0_1("v0.1"),

    /**
     * Differences between other versions:
     * - Spaces allowed in cells, for example in sample names with spaces in them ex: "C507/FG-CR 6".
     * - Fixed width fixed for floating point values
     */
    V0_2("v0.2"),

    /*
    * Differences between v0.x
    * - Added table and report headers
    * - Headers changed format, include the number of tables, rows, and metadata for gathering
    * - IS GATHERABLE
    */
    V1_0("v1.0"),

    /*
    * Differences between v1.0
    * - column numbers in header reflect the actual count of columns
    * - primary keys are never displayed
    */
    V1_1("v1.1");

    private final String versionString;

    private ReportVersion(String versionString) {
        this.versionString = versionString;
    }

    @Override
    public String toString() {
        return versionString;
    }

    public boolean equals(ReportVersion that) {
        return (versionString.equals(that.versionString));
    }

    /**
     * Returns the Report Version from the file header.
     *
     * @param header Header from the file starting with ##:Report.v[version]
     * @return The version as an enum.
     */
    public static ReportVersion fromHeader(String header) {
        if ( header == null )
            throw new UserException.BadInput("The report has no version specified in the header");

        if (header.startsWith("##:Report.v0.1 "))
            return ReportVersion.V0_1;

        if (header.startsWith("##:Report.v0.2 "))
            return ReportVersion.V0_2;

        if (header.startsWith("#:Report.v1.0"))
            return ReportVersion.V1_0;

        if (header.startsWith("#:Report.v1.1"))
            return ReportVersion.V1_1;

        throw new UserException.BadInput("The report has an unknown/unsupported version in the header: " + header);
    }
}
