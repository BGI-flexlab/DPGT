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
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.mapreduce.annotator;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

public class AnnotatorOptions extends GaeaOptions implements HadoopOptions{

    //the results file format:  tsv,vcf
    public enum OutputFormat {
        TSV, VCF
    }

    //the results file format:  tsv,vcf
    public enum RunStep {
        ANN, SORT, ALL
    }

    private final static String SOFTWARE_NAME = "Annotator";
    private final static String SOFTWARE_VERSION = "1.0";

    private String configFile; //用户配置文件
    private OutputFormat outputFormat;
    private String tmpPath;
    private String outputPath;
    private String inputFilePath;

    private String referenceSequencePath; //参考序列gaeaindex
    private RunStep runStep = RunStep.ALL;

    private boolean multiOutput = false;
    private boolean databaseCache = false;
    private boolean useDatabaseCache = false;
    private boolean verbose = false;
    private boolean debug = false;

    private int reducerNum;

    public AnnotatorOptions() {
        addOption("i", "input",      true,  "input file(VCF). [request]", true);
        addOption("o", "output",     true,  "output file of annotation results(.gz) [request]", true);
        addOption("c", "config",     true,  "config file. [request]", true);
        addOption("r", "reference",  true,  "indexed reference sequence file list [request]", true);
        addOption("O", "outputFormat", true,  "output file foramt (TSV, VCF) [TSV].");
        addOption("M", "MultiOutput", true,  "output to multi file when have more than one sample");
        addOption("d","doDatabaseCache",    false, "cache annotation database");
        addOption("u","useDatabaseCache",    false, "use cached annotation database");
        addOption("s","runStep",    true, "only run annotation step (ANN), only run sort step (SORT), run all steps (ALL) [ALL]");
        addOption(null,"verbose",    false, "display verbose information.");
        addOption(null,"debug",      false, "for debug.");
        addOption("h", "help",       false, "help information.");
        addOption("R", "reducer", true, "reducer numbers [30]");

        FormatHelpInfo(SOFTWARE_NAME,SOFTWARE_VERSION);
    }

    @Override
    public void parse(String[] args) {
        try {
            cmdLine = parser.parse(options, args);
            if(cmdLine.hasOption("h")) {
                helpInfo.printHelp("Options:", options, true);
                System.exit(1);
            }
        } catch (ParseException e) {
            helpInfo.printHelp("Options:", options, true);
            System.exit(1);
        }


        reducerNum = getOptionIntValue("R", 30);

        SimpleDateFormat df = new SimpleDateFormat("yyyyMMddHHmmss");
        tmpPath = "/user/" + System.getProperty("user.name") + "/annotmp-" + df.format(new Date());

        setInputFilePath(cmdLine.getOptionValue("input"));
        setConfigFile(cmdLine.getOptionValue("config"));
        setOutputFormat(OutputFormat.valueOf(cmdLine.getOptionValue("outputFormat", "TSV")));
        setReferenceSequencePath(cmdLine.getOptionValue("reference",""));
        setMultiOutput(getOptionBooleanValue("multiOutput", false));
        setOutputPath(cmdLine.getOptionValue("output"));
        setDatabaseCache(getOptionBooleanValue("doDatabaseCache", false));
        setUseDatabaseCache(getOptionBooleanValue("useDatabaseCache", false));
        setRunStep(getOptionValue("runStep", "ALL"));
        setVerbose(getOptionBooleanValue("verbose", false));
        setDebug(getOptionBooleanValue("debug", false));
    }

    @Override
    public void setHadoopConf(String[] args, Configuration conf) {
        conf.setStrings("args", args);
    }

    @Override
    public void getOptionsFromHadoopConf(Configuration conf) {
        String[] args = conf.getStrings("args");
        this.parse(args);
    }

    public int getReducerNum() {
        return reducerNum;
    }

    private String formatPath(String p) {
        if (p.startsWith("/")) {
            p = "file://" + p;
        }else if (p.startsWith(".")) {
            p = "file://" + new File(p).getAbsolutePath();
        }
        return p;
    }

    public String getOutputPath() {
        return outputPath;
    }

    public void setOutputPath(String outputPath) {
        this.outputPath = outputPath;
    }

    public String getTmpPath() {
        return tmpPath;
    }

    public void setTmpPath(String tmpPath) {
        this.tmpPath = tmpPath;
    }

    public String getInputFilePath() {
        return inputFilePath;
    }

    public void setInputFilePath(String inputFilePath) {
        this.inputFilePath = inputFilePath;
    }

    public String getConfigFile() {
        return configFile;
    }

    public void setConfigFile(String configFile) {
        this.configFile = formatPath(configFile);
    }

    public String getReferenceSequencePath() {
        return referenceSequencePath;
    }

    public void setReferenceSequencePath(String referenceSequencePath) {
        this.referenceSequencePath = formatPath(referenceSequencePath);
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public boolean isDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public String getOutput() {
        return formatPath(outputPath);
    }

    public String getInput() {
        return inputFilePath;
    }

    public boolean isMultiOutput() {
        return multiOutput;
    }

    public void setMultiOutput(boolean multiOutput) {
        this.multiOutput = multiOutput;
    }

    public OutputFormat getOutputFormat() {
        return outputFormat;
    }

    public void setOutputFormat(OutputFormat outputFormat) {
        this.outputFormat = outputFormat;
    }

    public boolean isDatabaseCache() {
        return databaseCache;
    }

    public void setDatabaseCache(boolean databaseCache) {
        this.databaseCache = databaseCache;
    }

    public boolean isUseDatabaseCache() {
        return useDatabaseCache;
    }

    public void setUseDatabaseCache(boolean useDatabaseCache) {
        this.useDatabaseCache = useDatabaseCache;
    }

    public RunStep getRunStep() {
        return runStep;
    }

    public void setRunStep(String runStep) {
        if(runStep.equalsIgnoreCase("ANN"))
            this.runStep = RunStep.ANN;
        else if(runStep.equalsIgnoreCase("SORT"))
            this.runStep = RunStep.SORT;
        else
            this.runStep = RunStep.ALL;
    }
}
