package org.bgi.flexlab.dpgt.jointcalling;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.Serializable;
import java.lang.Character;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.commons.cli.*;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.IntervalSetRule;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.bgi.flexlab.dpgt.utils.SimpleIntervalUtils;


public class JointCallingSparkOptions implements Serializable {
    private static final Logger logger = LoggerFactory.getLogger(JointCallingSparkOptions.class);
    // argument parser
    private Options options = new Options();
    private CommandLine cmdLine;
    private CommandLineParser parser = new PosixParser();
    private HelpFormatter helpFormatter = new HelpFormatter();

    // program meta data
    private static final String VERSION = "0.2.0.0";

    private ReferenceDataSource referenceDataSrc = null;
    private SAMSequenceDictionary sequenceDict = null;
    private List<SimpleInterval> targetIntervals = new ArrayList<>();

    // options
    public String input = null;             // input gvcf list in a file
    public String output = null;            // output directory
    public String reference = null;         // reference fasta file name
    public String[] targetRegionStrs = null;
    public String dbsnp = null;             // dbsnp vcf file

    public int jobs = 10;                   // number of spark exculators
    public int numCombinePartitions = 10;   // number of partitions for combining gvcfs.
    private static final String WINDOW = "300M";
    public int window = kmgStringToInt(WINDOW);   // window size for each combine-genotype cycle

    // genotype arguments
    public GenotypeCalculationArgumentCollection genotypeArguments = new GenotypeCalculationArgumentCollection();

    public boolean deleteIntermediateResults = true;

    public JointCallingSparkOptions() {
        addOption("i", "input", true, "input gvcf list in a file.", true, "FILE");
        addOption("o", "output", true, "output directory.", true, "DIR");
        addOption("r", "reference", true, "reference fasta file.", true, "FILE");
        addOption("l", "target-regions", true, "target regions to process.", false, "STRING");
        addOption("j", "jobs", true, "number of spark exculators.", false, "INT");
        addOption("n", "num-combine-partitions", true, "number of partitions for combining gvcfs. [10]", false, "INT");
        addOption("w", "window", true, "window size for each combine-genotype cycle. [300M]", false, "INT");
        addOption("d", "delete", true, "delete combine and genotype gvcf intermediate results. Possible values: true, false. [true]", false, "Boolean");
        addOption("s", "stand-call-conf", true, "the minimum phred-scaled confidence threshold at which variants should be called. [30.0]", false, "FLOAT");
        addOption(null, "dbsnp", true, "dbsnp vcf file for annotation.", false, "FILE");
        addOption(null, "use-old-qual-calculator", true, "use the old AF model. Possible values: true, false. [true]", false, "Boolean");
        addOption(null,"heterozygosity", true, "heterozygosity value used to compute prior likelihoods for any locus. [0.001]", false, "FLOAT");
        addOption(null,"indel-heterozygosity", true, "heterozygosity for indel calling. [1.25E-4]", false, "FLOAT");
        addOption(null,"heterozygosity-stdev", true, "standard deviation of heterozygosity for SNP and indel calling. [0.01]", false, "FLOAT");
        addOption(null,"ploidy", true, "ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy). [2]", false, "INT");
        addOption("h", "help", false, "print this message and exit.", false, null);
        addOption("V", "version", false, "show version.", false, null);
    }

    public void parse(String[] args) {
        formatHelpInfo();

        // check help and version
        for (String s : args) {
            if (s.equals("-h") || s.equals("--help")) {
                helpFormatter.printHelp("Options: ", options);
                System.exit(0);
            }
            if (s.equals("-V") || s.equals("--version")) {
                System.out.println("version: " + VERSION);
                System.exit(0);
            }
        }
        
        try {
            cmdLine = parser.parse(options, args);
        } catch (ParseException e) {
            logger.error("{}\n", e.getMessage());
            helpFormatter.printHelp("Options:", options);
            System.exit(1);
        }

        if (args.length == 0 || getOptionFlagValue("h")) {
            helpFormatter.printHelp("Options:", options);
            return;
        }

        if (getOptionFlagValue("V")) {
            System.out.println("DPGT Version: " + VERSION);
            return;
        }

        this.input = getOptionValue("i", null);
        this.output = getOptionValue("o", null);
        this.reference = getOptionValue("r", null);
        this.targetRegionStrs = getOptionValues("l", null);
        this.jobs = getOptionIntValue("j", this.jobs);
        this.numCombinePartitions = getOptionIntValue("n", this.numCombinePartitions);
        this.window = kmgStringToInt(getOptionValue("w", WINDOW));
        this.dbsnp = getOptionValue("dbsnp", this.dbsnp);
        this.deleteIntermediateResults = getOptionBooleanValue("d", this.deleteIntermediateResults);
        this.genotypeArguments.STANDARD_CONFIDENCE_FOR_CALLING = getOptionDoubleValue("s", this.genotypeArguments.STANDARD_CONFIDENCE_FOR_CALLING);
        this.genotypeArguments.useOldAFCalculator = getOptionBooleanValue("use-old-qual-calculator", true);  // note that we use old af calculator by default
        this.genotypeArguments.snpHeterozygosity = getOptionDoubleValue("heterozygosity", this.genotypeArguments.snpHeterozygosity);
        this.genotypeArguments.indelHeterozygosity = getOptionDoubleValue("indel-heterozygosity", this.genotypeArguments.snpHeterozygosity);
        this.genotypeArguments.heterozygosityStandardDeviation = getOptionDoubleValue("heterozygosity-stdev", this.genotypeArguments.heterozygosityStandardDeviation);
        this.genotypeArguments.samplePloidy = getOptionIntValue("ploidy", this.genotypeArguments.samplePloidy);

        this.referenceDataSrc = ReferenceDataSource.of(Paths.get(this.reference));
        this.sequenceDict = referenceDataSrc.getSequenceDictionary();
        if (this.targetRegionStrs != null) {
            GenomeLocSortedSet genomeLocs = IntervalUtils.loadIntervals(Arrays.asList(this.targetRegionStrs),
                IntervalSetRule.UNION, IntervalMergingRule.ALL, 0, new GenomeLocParser(this.sequenceDict));
            for (GenomeLoc loc: genomeLocs) {
                this.targetIntervals.add(new SimpleInterval(loc));
            }
        }
    }

    private void addOption(final String opt, final String longOpt, final boolean hasArg, final String description, final boolean isRequired, final String argName) {
        Option option = new Option(opt, longOpt, hasArg, description);
        option.setRequired(isRequired);
        if (hasArg) {
            option.setArgName(argName);
        }
        options.addOption(option);
    }

    public void formatHelpInfo() {
        helpFormatter.setSyntaxPrefix("DPGT: Distributed Population Genetics analysis Tools\nVersion: "+VERSION+"\n\n");
        helpFormatter.setWidth(HelpFormatter.DEFAULT_WIDTH*2);
        helpFormatter.setOptionComparator(null);
    }

    /**
     * convert string represent of size(eg 1k to 1000, 1m to 1000000, 1g to 1000000000) to int
     * @return
     */
    public static int kmgStringToInt(final String kmg) {
        final char unit = kmg.charAt(kmg.length() - 1);
        String numStr = null;
        int num = 0;
        if (Character.isDigit(unit)) {
            numStr = kmg.substring(0, kmg.length());
            try {
                num = Integer.parseInt(numStr);
            } catch (NumberFormatException e) {
                logger.error("Failed to convert string represent of size to int, {}", e.getMessage());
                System.exit(1);
            }
            return num;
        }

        numStr = kmg.substring(0, kmg.length() - 1);
        try {
            num = Integer.parseInt(numStr);
        } catch (NumberFormatException e) {
            logger.error("Failed to convert string represent of size to int, {}", e.getMessage());
            System.exit(1);
        }
        switch (unit) {
            case 'k':
            case 'K':
                num *= 1000;
                break;
            case 'm':
            case 'M':
                num *= 1000000;
                break;
            case 'g':
            case 'G':
                num *= 1000000000;
                break;
            default:
                logger.error("Invalid unit: {}", unit);
                System.exit(1);
        }

        return num;
    }

    public List<SimpleInterval> getIntervalsToTravers() {
        ArrayList<SimpleInterval> result = new ArrayList<>();
        if (this.targetIntervals == null) {
            for (SAMSequenceRecord loc: this.sequenceDict.getSequences()) {
                result.addAll(SimpleIntervalUtils.splitIntervalBySize(new SimpleInterval(loc), this.window));
            }
            return result;
        } else {
            for (SimpleInterval i: this.targetIntervals) {
                result.addAll(SimpleIntervalUtils.splitIntervalBySize(i, this.window));
            }
        }
        return result;
    }

    protected String[] getOptionValues(String opt, String[] defaultValue) {
        if (cmdLine.hasOption(opt)) 
            return cmdLine.getOptionValues(opt);
        return defaultValue;
    }

    protected String getOptionValue(String opt, String defaultValue) {
		if (cmdLine.hasOption(opt))
			return cmdLine.getOptionValue(opt);
		return defaultValue;
	}

	protected int getOptionIntValue(String opt, int defaultValue) {
		if (cmdLine.hasOption(opt))
			return Integer.parseInt(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

    protected boolean getOptionFlagValue(String opt) {
        if (cmdLine.hasOption(opt))
			return true;
		return false;
    }

	protected boolean getOptionBooleanValue(String opt, boolean defaultValue) {
		if (cmdLine.hasOption(opt)) {
            if (cmdLine.getOptionValue(opt).equals("true")) {
                return true;
            } else if (cmdLine.getOptionValue(opt).equals("false")) {
                return false;
            }
        }
        return defaultValue;
	}

	protected double getOptionDoubleValue(String opt, double defaultValue) {
		if (cmdLine.hasOption(opt))
			return Double.parseDouble(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected long getOptionLongValue(String opt, long defaultValue) {
		if (cmdLine.hasOption(opt))
			return Long.parseLong(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected byte getOptionByteValue(String opt, byte defaultValue) {
		if (cmdLine.hasOption(opt))
			return Byte.parseByte(cmdLine.getOptionValue(opt));
		return defaultValue;
	}

	protected short getOptionShortValue(String opt, short defaultValue) {
		if (cmdLine.hasOption(opt))
			return Short.parseShort(cmdLine.getOptionValue(opt));
		return defaultValue;
	}
}
