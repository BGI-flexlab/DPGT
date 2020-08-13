package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util;

import java.util.List;

import org.bgi.flexlab.gaea.util.StringUtils;

public final class AnnotationUtils {
    private AnnotationUtils(){}

    /**
     * Helper function to parse the list into the annotation string
     * @param valueList the ArrayList returned from StrandBiasBySample.annotate()
     * @return the array used by the per-sample Strand Bias annotation
     */
    public static String encodeValueList(final List<Double> valueList, final String precisionFormat ) {
        String[] outputList = new String[valueList.size()];
        int i =0;
        for (Double d : valueList) {
            outputList[i++] = String.format(precisionFormat, d);
        }
        return StringUtils.join(outputList, ",");
    }

    /**
     * Helper function to convert a List of Strings to a comma-separated String
     * @param stringList the ArrayList with String data
     * @return a comma-separated String
     */
    public static String encodeStringList( final List<String> stringList) {
    	String[] outputList = new String[stringList.size()];
        int i =0;
        for (String d : stringList) {
            outputList[i++] = d;
        }
        return StringUtils.join(outputList, ",");
    }

}
