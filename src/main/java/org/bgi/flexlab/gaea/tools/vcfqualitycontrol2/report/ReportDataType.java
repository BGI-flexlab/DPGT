package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.report;

import java.util.EnumSet;
import java.util.LinkedHashMap;
import java.util.Map;

public enum ReportDataType {
	 /**
     * The null type should not be used.
     */
    Null("Null"),

    /**
     * The default value when a format string is not present
     */
    Unknown("Unknown"),

    /**
     * Used for boolean values. Will display as true or false in the table.
     */
    Boolean("%[Bb]"),

    /**
     * Used for char values. Will display as a char so use printable values!
     */
    Character("%[Cc]"),

    /**
     * Used for float and double values. Will output a decimal with format %.8f unless otherwise specified.
     */
    Decimal("%.*[EeFf]"),

    /**
     * Used for int, byte, short, and long values. Will display the full number by default.
     */
    Integer("%[Dd]"),

    /**
     * Used for string values. Displays the string itself.
     */
    String("%[Ss]");

    private final java.lang.String dataTypeString;

    private ReportDataType(java.lang.String dataTypeString) {
        this.dataTypeString = dataTypeString;
    }

    private static final Map<java.lang.String, ReportDataType> lookup = new LinkedHashMap<>();

    static {
        for (ReportDataType s : EnumSet.allOf(ReportDataType.class))
            lookup.put(s.dataTypeString, s);
    }


    @Override
    public java.lang.String toString() {
        return this.dataTypeString;
    }

    /**
     * Returns a report data type from the Object specified. It looks through the list of acceptable classes and
     * returns the appropriate data type.
     *
     * @param object the object ot derive the data type from
     * @return the appropriate data type
     */
    public static ReportDataType fromObject(Object object) {
        ReportDataType value;
        if (object instanceof java.lang.Boolean) {
            value = ReportDataType.Boolean;

        } else if (object instanceof java.lang.Character) {
            value = ReportDataType.Character;

        } else if (object instanceof Float ||
                object instanceof Double) {
            value = ReportDataType.Decimal;

        } else if (object instanceof java.lang.Integer ||
                object instanceof Long ||
                object instanceof Short ||
                object instanceof Byte) {
            value = ReportDataType.Integer;

        } else if (object instanceof java.lang.String) {
            value = ReportDataType.String;

        } else {
            value = ReportDataType.Unknown;
        }
        return value;
    }

    /**
     * Returns a report data type from the format string specified. It uses regex matching from the enumerated
     * Strings.
     *
     * @param format the format string to derive the data type from
     * @return the appropriate data type
     */
    public static ReportDataType fromFormatString(java.lang.String format) {
        if (format.equals(""))
            return Unknown;
        for (ReportDataType type : lookup.values()) {
            if (format.matches(type.toString()) )
                return type;
        }
        return Unknown;
    }

    /**
     * Converts an input String to the appropriate type using the data type. Used for parsing loading a GATK report from
     * file.
     *
     * @param obj The input string
     * @return an object that matches the data type.
     */
    Object Parse(Object obj) {
        if (obj instanceof java.lang.String) {
            java.lang.String str = obj.toString();
            switch (this) {
                case Decimal:
                    return Double.parseDouble(str);
                case Boolean:
                    return java.lang.Boolean.parseBoolean(str);
                case Integer:
                    return Long.parseLong(str);
                case String:
                    return str;
                case Character:
                    return str.toCharArray()[0];
                default:
                    return str;
            }
        } else
            return null;
    }
}
