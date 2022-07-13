package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import java.nio.file.Paths;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

public class CombineVCFHeader {
    public Iterator<String> call(Iterator<String> vcfpathIter, final String output) {
        
        int inputIndex=0;
        LinkedHashSet<VCFHeader> headers=new LinkedHashSet<>();
        TreeSet<String> samples = new TreeSet<>();
        VCFHeader mergeHeader=null;
        while(vcfpathIter.hasNext()){
            String vcfpath=vcfpathIter.next();
            VCFFileReader reader =new VCFFileReader(Paths.get(vcfpath));
            inputIndex++;
            headers.add(reader.getHeader());
            samples.addAll(reader.getHeader().getSampleNamesInOrder());
            if(inputIndex>=1000){
                VCFHeader vcfHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true),samples);
                headers.clear();
                headers.add(vcfHeader);
            }
            reader.close();
        }

        mergeHeader=new VCFHeader(VCFUtils.smartMergeHeaders(headers, true), samples);
        headers.clear();

        VariantContextWriter writer = buildVCFHeaderWriter(output);
        writer.writeHeader(mergeHeader);
        writer.close();
        ArrayList<String> returnValue=new ArrayList<>();
        returnValue.add(output);
        return returnValue.iterator();
    }

    private VariantContextWriter buildVCFHeaderWriter(final String output) {
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(output);
        builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
        return builder.build();
    }
}
