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

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import java.util.ArrayList;
import java.util.List;
import java.lang.String;
import java.nio.file.Paths;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A MultiVariantSyncReader is a reader for reading one variant record from each
 * of the input vcf file and returning a list of variant records at a time.
 * The input vcfs are assumed to be vcfs for the same genome sites.
 */
public class MultiVariantSyncReader {
    private static final Logger logger = LoggerFactory.getLogger(MultiVariantSyncReader.class);
    private List<String> vcfpaths = null;
    private ArrayList<VCFFileReader> vcfReaders = null;
    private ArrayList<VariantContext> frontVCs = null;
    private ArrayList<CloseableIterator<VariantContext>> vcfIters = null;
    private VCFHeader header = null;
    private ArrayList<VariantContext> vcs = new ArrayList<>();

    public MultiVariantSyncReader() {}

    public void open(final List<String> vcfpaths, final String vcfHeader) {
        this.vcfpaths = vcfpaths;
        vcfReaders = new ArrayList<>();
        vcfIters = new ArrayList<>();
        for (String path: vcfpaths) {
            VCFFileReader reader = new VCFFileReader(Paths.get(path));
            vcfReaders.add(reader);
            vcfIters.add(reader.iterator());
        }
        VCFFileReader headerReader = new VCFFileReader(Paths.get(vcfHeader));
        header = headerReader.getHeader();
        headerReader.close();
        vcs = new ArrayList<>(vcfpaths.size());
    }

    public void close() {
        if (vcfIters != null) {
            for (CloseableIterator<VariantContext> iter: vcfIters) {
                iter.close();
            }
        }
        if (vcfReaders != null) {
            for (VCFFileReader reader: vcfReaders) {
                reader.close();
            }
        }
    }

    public void query(final String chrom, final int start, final int end) {
        int n = 0;
        for (VCFFileReader reader: vcfReaders) {
            vcfIters.set(n, reader.query(chrom, start, end));
            ++n;
        }
        frontVCs = getFrontVCs(start);
    }

    public void query(final Locatable locatable) {
        int n = 0;
        for (VCFFileReader reader: vcfReaders) {
            vcfIters.set(n, reader.query(locatable));
            ++n;
        }
        frontVCs = getFrontVCs(locatable.getStart());
    }

    public ArrayList<VariantContext> read() {
        if (frontVCs != null && !frontVCs.isEmpty()) {
            ArrayList<VariantContext> result = new ArrayList<>(frontVCs);
            frontVCs.clear();
            return result;
        }
        vcs.clear();
        CloseableIterator<VariantContext> firstIter = vcfIters.get(0);
        if (firstIter.hasNext()) {
            int n = 0;
            for (CloseableIterator<VariantContext> iter: vcfIters) {
                try {
                    vcs.add(iter.next());
                } catch (Exception e) {
                    logger.error("The vcf file {} has fewer lines than the first vcf file {}, {}",
                        vcfpaths.get(n), vcfpaths.get(0), e.getMessage());
                    System.exit(1);
                }
                ++n;
            }
        }
        return vcs;
    }

    public VCFHeader getVCFHeader() {
        return header;
    }

    /**
     * get variants at the front of query region, variant start outside of query
     * start (variant.start < query.start) will be skipped
     * @return
     */
    private ArrayList<VariantContext> getFrontVCs(int start) {
        ArrayList<VariantContext> results = new ArrayList<>();
        for (CloseableIterator<VariantContext> iter: vcfIters) {
            while(iter.hasNext()) {
                VariantContext vc = iter.next();
                if (vc.getStart() >= start) {
                    results.add(vc);
                    break;
                }
            }
        }
        return results;
    }

}
