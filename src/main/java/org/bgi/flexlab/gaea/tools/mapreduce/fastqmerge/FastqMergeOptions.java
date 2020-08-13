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
package org.bgi.flexlab.gaea.tools.mapreduce.fastqmerge;


import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;

public class FastqMergeOptions extends GaeaOptions implements HadoopOptions {

    private final static String SOFTWARE_NAME = "FastqMerge";
    private final static String SOFTWARE_VERSION = "1.0";

    private String fq1out;
    private String fq2out;
    private String input;

    private boolean verbose = false;
    private boolean debug = false;

    public FastqMergeOptions() {

        addOption("i", "input",      true,  "input dir. [request]", true);
        addOption("1", "fq1out",     true,  "read1 output file [request]", true);
        addOption("2", "fq2out",     true,  "read2 output file [null]");
        addOption(null,"verbose",    false, "display verbose information.");
        addOption(null,"debug",      false, "for debug.");
        addOption("h", "help",       false, "help information.");

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

        setInput(cmdLine.getOptionValue("input"));
        setFq1out(cmdLine.getOptionValue("fq1out"));
        setFq2out(cmdLine.getOptionValue("fq2out"));
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

    public String getInput() {
        return input;
    }

    public void setInput(String input) {
        this.input = input;
    }

    public String getFq1out() {
        return fq1out;
    }

    public void setFq1out(String fq1out) {
        this.fq1out = fq1out;
    }

    public String getFq2out() {
        return fq2out;
    }

    public void setFq2out(String fq2out) {
        this.fq2out = fq2out;
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

}
