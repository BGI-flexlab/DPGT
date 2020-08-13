package org.bgi.flexlab.gaea.tools.mapreduce.annotator.databaseload;

import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.mapreduce.options.HadoopOptions;
import org.bgi.flexlab.gaea.data.options.GaeaOptions;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RealignerExtendOptions;

public class DBNSFPToHbaseOptions extends GaeaOptions implements HadoopOptions {
    private final static String SOFTWARE_NAME = "DBNSFP2Hbase";
    private final static String SOFTWARE_VERSION = "1.0";

    // the path of hbase configuration core-site.xml hbase-site.xml
    private String hbaseConfig = null;

    // vcf file or directory path
    private String vcfInput = null;

    // hbase table name
    private String tableName = null;

    // hfile output path
    private String output = null;

    private int reduceNum = 1;

    public DBNSFPToHbaseOptions() {
        addOption("i", "input", true, "input directory", true);
        addOption("o", "output", true, "output directory", true);
        addOption("n", "reduceNumber", true, "reduce number");
        addOption("c", "config", true, "hbase configuration path", true);
        addOption("t", "table", true, "hbase table name", true);
        FormatHelpInfo(SOFTWARE_NAME, SOFTWARE_VERSION);
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

    @Override
    public void parse(String[] args) {
        try {
            cmdLine = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            FormatHelpInfo(RealignerExtendOptions.SOFTWARE_NAME, RealignerExtendOptions.SOFTWARE_VERSION);
            System.exit(1);
        }

        vcfInput = getOptionValue("i",null);
        output = getOptionValue("o",null);
        tableName = getOptionValue("t",null);
        hbaseConfig = getOptionValue("c",null);
        reduceNum = getOptionIntValue("n",1);

        checkArgs();
    }

    public void checkArgs(){
        if(vcfInput == null)
            throw new UserException.BadArgumentValueException("i","null");
        if(output == null)
            throw new UserException.BadArgumentValueException("o","null");
        if(hbaseConfig == null)
            throw new UserException.BadArgumentValueException("c","null");
        if(tableName == null)
            throw new UserException.BadArgumentValueException("t","null");
    }

    public String getTableName(){
        return this.tableName;
    }

    public String getInput(){
        return this.vcfInput;
    }

    public int getReducerNumber(){
        return this.reduceNum;
    }

    public String getHeaderOutput(){
        if(output.endsWith("/"))
            return this.output+"header";
        return this.output+"/header";
    }

    public String getHFileOutput(){
        if(output.endsWith("/"))
            return this.output+"hfile";
        return this.output+"/hfile";
    }

    public String getConfig(){
        if(this.hbaseConfig.endsWith("/"))
            return this.hbaseConfig;
        return this.hbaseConfig+"/";
    }
}
