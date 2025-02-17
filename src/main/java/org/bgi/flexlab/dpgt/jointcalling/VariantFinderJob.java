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
import java.util.List;
import java.util.BitSet;
import java.util.Arrays;
import org.bgi.flexlab.dpgt.utils.DPGTJob;
import org.bgi.flexlab.dpgt.utils.VariantSiteSetUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;


public class VariantFinderJob extends DPGTJob<BitSet> {
    private JavaRDD<String> vcfpathsRDDPartitionByJobs;
    private JointCallingSparkOptions jcOptions;
    private JavaSparkContext sc;
    private int idx;
    private SimpleInterval interval;
    public final File variantSiteDir;

    public VariantFinderJob(final JavaRDD<String> vcfpathsRDDPartitionByJobs, final JointCallingSparkOptions jcOptions,
        final JavaSparkContext sc, final int idx, final SimpleInterval interval) {
        this.vcfpathsRDDPartitionByJobs = vcfpathsRDDPartitionByJobs;
        this.jcOptions = jcOptions;
        this.sc = sc;
        this.idx = idx;
        this.interval = interval;
        this.variantSiteDir = new File(this.jcOptions.output+"/"+JointCallingSparkConsts.VARIANT_SITE_PREFIX+this.idx);
        if (!this.variantSiteDir.exists()) {
            this.variantSiteDir.mkdirs();
        }
        this.stateFile = this.jcOptions.output + "/" + JointCallingSparkConsts.JOB_STATE + "/" + JointCallingSparkConsts.VARIANT_SITE_PREFIX + this.idx + ".json";
    }


    public BitSet work() {
        final String variantSitePrefix = variantSiteDir.getAbsolutePath()+"/"+JointCallingSparkConsts.VARIANT_SITE_PREFIX;
        List<String> variantSiteFiles = vcfpathsRDDPartitionByJobs
            .mapPartitionsWithIndex(new VariantSiteFinderSparkFunc(variantSitePrefix, interval.getContig(), interval.getStart()-1, interval.getEnd()-1), false)
            .filter(x -> {return !x.equals("null");})
            .collect();
        
        JavaRDD<String> variantSiteFilesRDD = sc.parallelize(variantSiteFiles, 1);
        List<BitSet> variantSiteSetDatas = variantSiteFilesRDD
            .mapPartitionsWithIndex(new CombineVariantSiteSetSparkFunc(variantSitePrefix), false)
            .collect();

        BitSet variantSiteSetData = null;
        if (!variantSiteSetDatas.isEmpty()) variantSiteSetData = variantSiteSetDatas.get(0);

        final String variantSiteFile = variantSiteDir.getAbsolutePath()+"/"+JointCallingSparkConsts.VARIANT_SITE_PREFIX+JointCallingSparkConsts.VARIANT_SITE_SUFFIX.substring(1);
        VariantSiteSetUtils.writeVariantSiteSet(variantSiteSetData, variantSiteFile);
        this.jobState.outPutFiles.put(0, Arrays.asList(new String[]{variantSiteFile}));

        // delete intermediate variant site files
        for (String f: variantSiteFiles) {
            File file = new File(f);
            file.delete();
        }

        return variantSiteSetData;
    }

    public BitSet load() {
        return VariantSiteSetUtils.loadVariantSiteSet(this.jobState.outPutFiles.get(0).get(0));
    }
}
