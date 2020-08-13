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
package org.bgi.flexlab.gaea.tools.mapreduce.vcfstats;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;

import java.text.SimpleDateFormat;
import java.util.Date;

public class VCFStatsOptions extends GaeaOptions implements HadoopOptions {

    private final static String SOFTWARE_NAME = "VCFStats";
    private final static String SOFTWARE_VERSION = "1.0";

    private String tmpPath;
    private String outputPath;
    private String input;
    private String dbsnpFile; //vcf.gz, must be indexed

    private String referenceSequencePath; //参考序列gaeaindex

    private boolean countVarLength = false;
    private boolean verbose = false;
    private boolean debug = false;

    private int reducerNum;

    public VCFStatsOptions() {
        addOption("i", "input",      true,  "input file(VCF). [request]", true);
        addOption("o", "output",     true,  "output file [request]", true);
        addOption("d", "dbsnp",     true,  "dbsnp file(.vcf.gz), must be indexed [null]");
        addOption("r", "reference",  true,  "indexed reference sequence file list [request]", true);
        addOption("c", "countVarLength",  false,  "count variant length");
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
        tmpPath = "/user/" + System.getProperty("user.name") + "/vcfsummarytmp-" + df.format(new Date());

        setInput(cmdLine.getOptionValue("input"));
        setDbsnpFile(cmdLine.getOptionValue("dbsnp"));
        setReferenceSequencePath(cmdLine.getOptionValue("reference",""));
        setOutputPath(cmdLine.getOptionValue("output"));
        setCountVarLength(getOptionBooleanValue("c", false));
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

    public String getInput() {
        return input;
    }

    public void setInput(String input) {
        this.input = input;
    }

    public String getReferenceSequencePath() {
        return referenceSequencePath;
    }

    public void setReferenceSequencePath(String referenceSequencePath) {
        this.referenceSequencePath = referenceSequencePath;
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public boolean isCountVarLength() {
        return countVarLength;
    }

    public void setCountVarLength(boolean countVarLength) {
        this.countVarLength = countVarLength;
    }

    public boolean isDebug() {
        return debug;
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    public String getDbsnpFile() {
        return dbsnpFile;
    }

    public void setDbsnpFile(String dbsnpFile) {
        this.dbsnpFile = dbsnpFile;
    }
}
