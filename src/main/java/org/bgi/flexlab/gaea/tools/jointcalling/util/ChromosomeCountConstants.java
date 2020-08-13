package org.bgi.flexlab.gaea.tools.jointcalling.util;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class ChromosomeCountConstants {

    public static final String[] keyNames = { VCFConstants.ALLELE_NUMBER_KEY, VCFConstants.ALLELE_COUNT_KEY, VCFConstants.ALLELE_FREQUENCY_KEY };

    public static final VCFInfoHeaderLine[] descriptions = {
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY),
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY),
            VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY) };
}
