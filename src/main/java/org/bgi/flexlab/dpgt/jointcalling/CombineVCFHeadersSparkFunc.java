package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.apache.spark.api.java.function.Function2;
import java.nio.file.Paths;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

public class CombineVCFHeadersSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String outdir;

    public CombineVCFHeadersSparkFunc(final String outdir) {
        this.outdir = outdir;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) throws Exception {
        ArrayList<String> vcfpaths = new ArrayList<>();
        while(vcfpathIter.hasNext()) {
            vcfpaths.add(vcfpathIter.next());
        }
        if (vcfpaths.isEmpty()) {
            ArrayList<String> returnValue=new ArrayList<>();
            returnValue.add("null");
            return returnValue.iterator();
        }
        String[] vcfpathsArray  = new String[vcfpaths.size()];
        vcfpaths.toArray(vcfpathsArray);

        String output = this.outdir + "/header." + idx + ".vcf.gz";
        VCFHeaderCombiner combiner = new VCFHeaderCombiner();
        combiner.Combine(vcfpathsArray, output);
        ArrayList<String> returnValue=new ArrayList<>();
        returnValue.add(output);
        return returnValue.iterator();
    }
}
