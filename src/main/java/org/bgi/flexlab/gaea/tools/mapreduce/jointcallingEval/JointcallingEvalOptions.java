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
package org.bgi.flexlab.gaea.tools.mapreduce.jointcallingEval;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Created by huangzhibo on 2017/12/12.
 */
public class JointcallingEvalOptions extends GaeaOptions implements HadoopOptions {
    private final static String SOFTWARE_NAME = "JointcallingEval";
    private final static String SOFTWARE_VERSION = "1.0";

    public enum Mode {
        SNP, INDEL, BOTH
    }

    private String tmpPath = null; //输出格式 txt,vcf
    private String outputPath = null;
    private String inputFilePath = null;
    private String baselineFile = null;
    private int reducerNum;
    private Mode mode;
    private boolean outputdiff;

    public JointcallingEvalOptions() {
        addOption("i", "input",      true,  "input file(VCF). [request]", true);
        addOption("b", "baseline",      true,  "baseline file(VCF). [request]", true);
        addOption("o", "output",     true,  "output file of annotation results(.gz) [request]", true);
        addOption("d", "outputdiff",     false,  "output diff variations (.gz) [false]");
        addOption("m", "mode",     true,  "only count SNP， INDEL or BOTH [BOTH]");
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
        setOutputPath(cmdLine.getOptionValue("output"));
        setBaselineFile(cmdLine.getOptionValue("baseline"));
        setMode(getOptionValue("m", "BOTH"));
        setOutputdiff(getOptionBooleanValue("d", false));
    }

    @Override
    public void setHadoopConf(String[] args, Configuration conf) {
        conf.setStrings("args", args);
        conf.set("inputFilePath", getInputFilePath());
        conf.set("outputPath", getOutputTmpPath());
    }

    @Override
    public void getOptionsFromHadoopConf(Configuration conf) {
        String[] args = conf.getStrings("args");
        this.parse(args);
    }

    private String formatPath(String p) {
        if (p.startsWith("/")) {
            p = "file://" + p;
        }else if (p.startsWith(".")) {
            p = "file://" + new File(p).getAbsolutePath();
        }
        return p;
    }

    public int getReducerNum() {
        return reducerNum;
    }

    public String getOutputTmpPath() {
        return formatPath(outputPath+"/tmp");
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

    public String getInput() {
        return inputFilePath;
    }

    public void setInputFilePath(String inputFilePath) {
        this.inputFilePath = inputFilePath;
    }

    public String getBaselineFile() {
        return baselineFile;
    }

    public void setBaselineFile(String baselineFile) {
        this.baselineFile = baselineFile;
    }

    public boolean isOutputdiff() {
        return outputdiff;
    }

    public void setOutputdiff(boolean outputdiff) {
        this.outputdiff = outputdiff;
    }

    public Mode getMode() {
        return mode;
    }

    public void setMode(String mode) {
        if(mode.equalsIgnoreCase("SNP"))
            this.mode = Mode.SNP;
        else if(mode.equalsIgnoreCase("INDEL"))
            this.mode = Mode.INDEL;
        else
            this.mode = Mode.BOTH;
    }
}
