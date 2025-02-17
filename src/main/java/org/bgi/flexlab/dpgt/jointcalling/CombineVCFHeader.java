/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
package org.bgi.flexlab.dpgt.jointcalling;

import java.util.*;
import java.util.Iterator;
import java.io.IOException;
import java.nio.file.Paths;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class CombineVCFHeader {
    private static final Logger logger = LoggerFactory.getLogger(CombineVCFHeader.class);
    public Iterator<String> call(Iterator<String> vcfpathIter, final String output) {
        
        int inputIndex=0;
        LinkedHashSet<VCFHeader> headers=new LinkedHashSet<>();
        TreeSet<String> samples = new TreeSet<>();
        VCFHeader mergeHeader=null;
        while(vcfpathIter.hasNext()){
            String vcfpath=vcfpathIter.next();
            VCFIterator reader = null;
            try {
                reader =new VCFIteratorBuilder().open(vcfpath);
            } catch (IOException e) {
                logger.error("{}", e.getMessage());
                System.exit(1);
            }
            inputIndex++;
            headers.add(reader.getHeader());
            samples.addAll(reader.getHeader().getSampleNamesInOrder());
            if(inputIndex>=1000){
                inputIndex = 0;
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
