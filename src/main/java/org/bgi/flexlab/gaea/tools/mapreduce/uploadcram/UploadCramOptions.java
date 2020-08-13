package org.bgi.flexlab.gaea.tools.mapreduce.uploadcram;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.util.GenericOptionsParser;
import org.apache.hadoop.util.LineReader;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.seqdoop.hadoop_bam.SAMFormat;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by huangzhibo on 2018/7/23.
 */
public class UploadCramOptions extends GaeaOptions implements HadoopOptions {

    private final static String SOFTWARE_NAME = "UploadCram";
    private final static String SOFTWARE_VERSION = "1.0";

    public UploadCramOptions() {
        // TODO Auto-generated constructor stub
        addOption("i", "input", true, "input file or dir", true);
        addOption("o", "output", true, "output dir", true);
        addOption("r", "ref", true, "reference file for read cram", true);
        addOption("h", "help", false, "help information");
    }

    private List<Path> inputs = new ArrayList<>();

    private String outputPath;

    private String localReferenceSequencePath;

    /**
     * regions in bed file
     */
    private String bedfile;

    @Override
    public void setHadoopConf(String[] args, Configuration conf) {
        // TODO Auto-generated method stub
        String[] otherArgs;
        try {
            otherArgs = new GenericOptionsParser(args).getRemainingArgs();
            conf.setStrings("args", otherArgs);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    @Override
    public void getOptionsFromHadoopConf(Configuration conf) {
        // TODO Auto-generated method stub
        String[] args = conf.getStrings("args");
        this.parse(args);
    }

    @Override
    public void parse(String[] args) {
        // TODO Auto-generated method stub
        try {
            cmdLine = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            printHelpInfotmation(SOFTWARE_NAME);
            System.exit(1);
        }

        try {
            parseInput( getOptionValue("i",null));
        } catch (IOException e) {
            throw new UserException(e.toString());
        }

        localReferenceSequencePath = getOptionValue("r", null);
        outputPath = getOptionValue("o", null);
        bedfile = getOptionValue("B", null);
    }


    private void parseInput(String input) throws IOException {
        Path path = new Path(input);
        FileSystem inFS = path.getFileSystem(new Configuration());
        if(inFS.isDirectory(path)){
            FileStatus[] fileStatuses = inFS.globStatus(new Path(input +"/part*"));
            for (FileStatus f: fileStatuses)
                inputs.add(f.getPath());
            return;
        }

        if(input.endsWith(".cram")){
            inputs.add(path);
        }
        else {
            LineReader reader = new LineReader(inFS.open(path));
            Text line = new Text();

            while(reader.readLine(line) > 0 && line.getLength() != 0) {
                inputs.add(new Path(line.toString()));
            }
            reader.close();
        }
    }

    public List<Path> getInputs() {
        return inputs;
    }

    public void setInputs(List<Path> inputs) {
        this.inputs = inputs;
    }

    /**
     * @return the outputPath
     */
    public String getOutputPath() {
        return outputPath;
    }

    /**
     * @param outputPath the outputPath to set
     */
    public void setOutputPath(String outputPath) {
        this.outputPath = outputPath;
    }


    public String getLocalReferenceSequencePath() {
        return localReferenceSequencePath;
    }

}
