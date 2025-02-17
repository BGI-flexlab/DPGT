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
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.bgi.flexlab.dpgt.utils.NativeLibraryLoader;


public class CombineGVCFsOnSitesSparkFunc implements Function2<Integer, Iterator<String>, Iterator<String>> {
    public String refpath;
    public String prefix;
    public SimpleInterval interval;
    Broadcast<byte[]> variantSiteSetBytesBc;

    static {
        NativeLibraryLoader.load();
    }

    public CombineGVCFsOnSitesSparkFunc(final String refpath,
        final String prefix, final SimpleInterval interval, Broadcast<byte[]> variantSiteSetBytesBc)
    {
        this.refpath = refpath;
        this.prefix = prefix;
        this.interval = interval;
        this.variantSiteSetBytesBc = variantSiteSetBytesBc;
    }

    @Override public Iterator<String> call(Integer idx, Iterator<String> vcfpathIter) {
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
        CombineGVCFsOnSites combiner = new CombineGVCFsOnSites();
        String outpath = prefix + idx + ".vcf.gz";
        // call native c++ function to combine gvcfs
        combiner.Combine(vcfpathsArray, refpath, outpath, variantSiteSetBytesBc.value(), interval.getContig(), interval.getStart()-1, interval.getEnd()-1);
        ArrayList<String> result=new ArrayList<>();
        result.add(outpath);
        return result.iterator();
    }
}
