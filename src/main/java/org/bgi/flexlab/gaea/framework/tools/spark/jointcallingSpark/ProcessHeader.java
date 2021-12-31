package org.bgi.flexlab.gaea.framework.tools.spark.jointcallingSpark;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.Path;
import org.apache.spark.api.java.function.Function2;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalWriter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class ProcessHeader implements Function2<Integer,Iterator<String>,Iterator<String>> {
    public JointCallingSparkOptions options;
    public ProcessHeader(JointCallingSparkOptions options){
        this.options=options;
    }
    @Override public Iterator<String> call(Integer v1, Iterator<String> v2) throws Exception {
        Logger log= LoggerFactory.getLogger(ProcessHeader.class);
        BufferedWriter vcfHeaderWriter=new BufferedWriter(new FileWriter(options.getOutDir()+"/headers/vcfPathName"+v1));
        int inputIndex=0;
        LinkedHashSet<VCFHeader> headers=new LinkedHashSet<>();
        TreeSet<String> sampleNames=new TreeSet<>();
        int timeIter=0;
        VCFHeader mergeHeader=null;
        while(v2.hasNext()){
            String a=v2.next();
            if(a.startsWith("file://")){
                a=a.substring(7);
            }
            File aFile=new File(a);
            if(timeIter%100==0){
                log.info("current sample index:\t"+timeIter);
            }
            timeIter++;
            VCFLocalLoader vcfLL=new VCFLocalLoader(a);
            vcfHeaderWriter.write(a+"\t"+aFile.getName()+"\t"+vcfLL.getHeader().getSampleNamesInOrder().get(0)+"\n");
            inputIndex++;
            headers.add(vcfLL.getHeader());
            sampleNames.addAll(vcfLL.getHeader().getSampleNamesInOrder());
            if(inputIndex>=1000){
                VCFHeader vcfHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),sampleNames);
                headers.clear();
                headers.add(vcfHeader);
            }
        }


        mergeHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),sampleNames);
        headers.clear();

        vcfHeaderWriter.close();
        VCFLocalWriter mergedVCFHeaderWriter=new VCFLocalWriter(options.getOutDir()+"/headers/vcfheader"+v1,false,false);
        mergedVCFHeaderWriter.writeHeader(mergeHeader);
        mergedVCFHeaderWriter.close();
        ArrayList<String> returnValue=new ArrayList<>();
        returnValue.add("done"+v1);
        return returnValue.iterator();
    }
}
