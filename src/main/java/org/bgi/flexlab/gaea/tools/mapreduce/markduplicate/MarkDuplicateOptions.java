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
package org.bgi.flexlab.gaea.tools.mapreduce.markduplicate;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.fs.PathFilter;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by huangzhibo on 2017/4/14.
 */
public class MarkDuplicateOptions extends GaeaOptions implements HadoopOptions {
    private final static String SOFTWARE_NAME = "MarkDuplicate";
    private final static String SOFTWARE_VERSION = "1.0";

    private String input;
    private ArrayList<Path> inputFileList;
    private int inputFormat;
    private String output;
    private int outputFormat;
    private boolean outputDupRead;
    private boolean isSE;
    private boolean removeSecond;
    private int reducerNum;
    private int windowSize;
    private int extendSize;
    FileSystem fs;

    public MarkDuplicateOptions() {
        addOption("i", "input", true, "input directory [required]", true);
        addOption("I", "inputFormat", true, "input Format. 0:BAM; 1:SAM [1]");
        addOption("o", "output", true, "output directory [required]", true);
        addOption("O", "outputFormat", true, "output Format. 0:BAM; 1:SAM [0]");
        addOption("D", "outputDupRead", false, "output Duplicates reads [true]");
        addOption("M", "removeSecond", false, "remove not primary and supplementary alignment reads [false]");
        addOption("S", "isSE", false, "input is SE data [false]");
        addOption("R", "reducer", true, "reducer numbers [30]");
        addOption("W", "windowSize", true, "window size that sharding the data [100000]");
        addOption("E", "extendSize", true, "The extend size (must greater than read length) [100]");
        addOption("h", "help", false, "print help information.");
        FormatHelpInfo(SOFTWARE_NAME,SOFTWARE_VERSION);

        inputFileList = new ArrayList<>();
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

        if(args.length == 0 || getOptionBooleanValue("h", false)) {
            printHelpInfotmation(SOFTWARE_NAME);
            System.exit(1);
        }

        input = getOptionValue("i", null);
        inputFormat = getOptionIntValue("I", 1);
        output = getOptionValue("o", null);
        outputFormat = getOptionIntValue("O", 0);
        outputDupRead = getOptionBooleanValue("D", true);
        removeSecond = getOptionBooleanValue("M", false);
        isSE = getOptionBooleanValue("S", false);
        reducerNum = getOptionIntValue("R", 30);
        windowSize = getOptionIntValue("W", 100000);
        extendSize = getOptionIntValue("E", 100);
    }

    @Override
    public void setHadoopConf(String[] args, Configuration conf) {
        conf.setStrings("args", args);
        Path p = new Path(this.getInput());
        try {
            fs = p.getFileSystem(conf);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        traversalInputPath(this.getInput());
    }

    @Override
    public void getOptionsFromHadoopConf(Configuration conf) {
        String[] args = conf.getStrings("args");
        this.parse(args);
    }

    private void traversalInputPath(String input)
    {
        Path path = new Path(input);
        PathFilter filter = file -> !file.getName().startsWith("_");
        try {
            if (!fs.exists(path)) {
                System.err.println("Input File Path is not exist! Please check -i var.");
                System.exit(-1);
            }
            if (fs.isFile(path)) {
                inputFileList.add(path);
            }else {
                FileStatus[] stats=fs.listStatus(path, filter);
                if(stats.length <= 0){
                    System.err.println("Input File Path is empty! Please check input : " +path.toString());
                    System.exit(-1);
                }

                for (FileStatus file : stats) {
                    Path filePath=file.getPath();

                    if (!fs.isFile(filePath)) {
                        String childPath=filePath.toString();
                        traversalInputPath(childPath);
                    }else {
                        inputFileList.add(filePath);
                    }
                }
            }
        }catch (IOException ioe) {
            throw new RuntimeException(ioe);
        }
    }

    public ArrayList<Path> getInputFileList(){
        return inputFileList;
    }

    public String getInput() {
        return input;
    }

    public SAMFormat getInputFormat(){
        return inputFormat == 0 ? SAMFormat.BAM :SAMFormat.SAM;
    }

    public String getOutput() {
        return output;
    }

    public int getOutputFormat() {
        return outputFormat;
    }

    public boolean isOutputDupRead() {
        return outputDupRead;
    }

    public boolean isSE() {
        return isSE;
    }

    public boolean isRemoveSecond() {
        return removeSecond;
    }

    public void setRemoveSecond(boolean removeSecond) {
        this.removeSecond = removeSecond;
    }

    public int getReducerNum() {
        return reducerNum;
    }

    public int getWindowSize() {
        return windowSize;
    }

    public int getExtendSize() {
        return extendSize;
    }
}
