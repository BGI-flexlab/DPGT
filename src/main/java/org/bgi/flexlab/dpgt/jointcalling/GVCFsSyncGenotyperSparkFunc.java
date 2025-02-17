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

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.spark.api.java.function.Function2;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;


public class GVCFsSyncGenotyperSparkFunc implements Function2<Integer, Iterator<SimpleInterval>, Iterator<String>> {
    public String refpath;
    public List<String> vcfpaths;
    public String vcfHeaderPath;
    public String prefix;
    public String dbsnpPath;
    public GenotypeCalculationArgumentCollection genotypeArgs;
    
    public GVCFsSyncGenotyperSparkFunc(final String refpath, final List<String> vcfpaths, final String vcfHeaderPath,
        final String prefix, final String dbsnpPath, final GenotypeCalculationArgumentCollection genotypeArgs) {
        this.refpath = refpath;
        this.vcfpaths = vcfpaths;
        this.vcfHeaderPath = vcfHeaderPath;
        this.prefix = prefix;
        this.dbsnpPath = dbsnpPath;
        this.genotypeArgs = genotypeArgs;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<SimpleInterval> intervalIter) {
        if (intervalIter.hasNext()) {
            SimpleInterval interval = intervalIter.next();
            String outpath = prefix + idx + ".vcf.gz";
            GVCFsSyncGenotyper genotyper = new GVCFsSyncGenotyper(refpath, vcfpaths, vcfHeaderPath, interval, outpath, dbsnpPath, genotypeArgs);
            genotyper.run();
            genotyper.stop();
            ArrayList<String> result=new ArrayList<>();
            result.add(outpath);
            return result.iterator();
        }
        ArrayList<String> result=new ArrayList<>();
        result.add("null");
        return result.iterator();
    }
}
