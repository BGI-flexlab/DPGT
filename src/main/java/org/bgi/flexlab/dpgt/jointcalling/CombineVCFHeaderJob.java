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

import java.io.File;
import java.util.Arrays;
import java.util.List;
import org.bgi.flexlab.dpgt.utils.DPGTJob;
import org.apache.spark.api.java.JavaRDD;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;


public class CombineVCFHeaderJob extends DPGTJob<String> {
    private JavaRDD<String> vcfpathsRDD;
    private JointCallingSparkOptions jcOptions;
    public File headerDir;

    public CombineVCFHeaderJob(final JavaRDD<String> vcfpathsRDD, final JointCallingSparkOptions jcOptions) {
        this.vcfpathsRDD = vcfpathsRDD;
        this.jcOptions = jcOptions;
        this.headerDir = new File(jcOptions.output + "/" + JointCallingSparkConsts.HEADER_DIR);
        if (!headerDir.exists()) {
            headerDir.mkdirs();
        }
        File jobStateDir = new File(jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE);
        if (!jobStateDir.exists()) {
            jobStateDir.mkdirs();
        }
        this.stateFile = jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.COMBINE_HEADER_STATE_FILE;
    }

    public String work() {
        List<String> combinedHeaders = vcfpathsRDD.mapPartitionsWithIndex(new CombineVCFHeadersSparkFunc(headerDir.getAbsolutePath()), false)
            .filter(x -> {return !x.equals("null");})
            .collect();

        // combine headers of each partition to generate combined header for all input vcfs
        VCFHeaderCombiner combiner = new VCFHeaderCombiner();
        String combinedHeader = headerDir.getAbsolutePath()+"/"+JointCallingSparkConsts.COMBINE_HEADER;
        String[] combinedHeadersArray = new String[combinedHeaders.size()];
        combinedHeaders.toArray(combinedHeadersArray);
        combiner.Combine(combinedHeadersArray, combinedHeader);

        // make combined header for genotyping by add genotyping specific headers
        GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper();
        VCFHeader combinedHeaderForGenotyping = genotyper.makeCombinedHeaderForGenotyping(combinedHeader, jcOptions.genotypeArguments);
        String combinedHeaderForGenotypingPath = headerDir.getAbsolutePath()+"/"+JointCallingSparkConsts.GENOTYPE_HEADER;
        VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(combinedHeaderForGenotypingPath);
        builder.setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF);
        VariantContextWriter writer = builder.build();
        writer.writeHeader(combinedHeaderForGenotyping);
        writer.close();

        this.jobState.outPutFiles.put(0, Arrays.asList(new String[]{combinedHeaderForGenotypingPath}));
        return combinedHeaderForGenotypingPath;
    }

    public String load() {
        return jobState.outPutFiles.get(0).get(0);
    }
}
