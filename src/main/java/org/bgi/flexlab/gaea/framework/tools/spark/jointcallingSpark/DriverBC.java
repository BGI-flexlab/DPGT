package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import scala.Serializable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

public class DriverBC implements Serializable {
    public DriverBC(String outputDir, ArrayList<ArrayList<String>> multiMapSampleNames, LinkedHashMap<String, Integer> chrIndex, JointCallingSparkOptions options, String vVcfPath, LinkedHashMap<String,Integer> sampleIndex, Map<String, String> pathSample, HashMap<Integer,Long> accumulateLength,
                    VCFHeader virtualHeader,VCFHeaderVersion version){
        this.outputDir=outputDir;
        this.multiMapSampleNames.addAll(multiMapSampleNames);
        this.chrIndex.putAll(chrIndex);
        this.options=options;
        this.vVcfPath=vVcfPath;
        this.sampleIndex.putAll(sampleIndex);
        this.pathSample.putAll(pathSample);
        this.accumulateLength.putAll(accumulateLength);
        this.virtualHeader=virtualHeader;
        this.version=version;
    }
    public String outputDir;
    public final static String INPUT_ORDER = "input.name.order";
    public final static String Window_File = "window.file.path";
    public ArrayList<ArrayList<String>> multiMapSampleNames=new ArrayList<>();
    public LinkedHashMap<String, Integer> chrIndex=new LinkedHashMap<>();
    public JointCallingSparkOptions options;

    static {
        new JointCallingSparkOptions();
    }

    public String vVcfPath;
    public LinkedHashMap<String,Integer> sampleIndex=new LinkedHashMap<>();
    public Map<String, String> pathSample = new HashMap<>();
    public HashMap<Integer,Long> accumulateLength=new HashMap<>();
    public VCFHeader virtualHeader;
    public VCFHeaderVersion version;
}
