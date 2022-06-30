package org.bgi.flexlab.dpgt;

import htsjdk.variant.vcf.VCFFileReader;
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
    private ArrayList<CloseableIterator<VariantContext>> vcfIters = null;

    public MultiVariantSyncReader() {}

    public void open(final List<String> vcfpaths) {
        this.vcfpaths = vcfpaths;
        vcfReaders = new ArrayList<>();
        vcfIters = new ArrayList<>();
        for (String path: vcfpaths) {
            VCFFileReader reader = new VCFFileReader(Paths.get(path));
            vcfReaders.add(reader);
            vcfIters.add(reader.iterator());
        }
    }

    public void close() {
        if (vcfIters != null) {
            for (CloseableIterator<VariantContext> iter: vcfIters) {
                iter.close();
            }
        } else if (vcfReaders != null) {
            for (VCFFileReader reader: vcfReaders) {
                reader.close();
            }
        }
    }

    public void query(final String chrom, final int start, final int end) {
        int n = 0;
        for (VCFFileReader reader: vcfReaders) {
            reader.query(chrom, start, end);
            vcfIters.set(n, reader.iterator());
            ++n;
        }
    }

    public void query(final Locatable locatable) {
        int n = 0;
        for (VCFFileReader reader: vcfReaders) {
            reader.query(locatable);
            vcfIters.set(n, reader.iterator());
            ++n;
        }
    }

    public ArrayList<VariantContext> read() {
        ArrayList<VariantContext> result = new ArrayList<>();
        CloseableIterator<VariantContext> firstIter = vcfIters.get(0);
        if (firstIter.hasNext()) {
            int n = 0;
            for (CloseableIterator<VariantContext> iter: vcfIters) {
                try {
                    result.add(iter.next());
                } catch (Exception e) {
                    logger.error("The vcf file %s has fewer lines than the first vcf file %s, %s",
                        vcfpaths.get(n), vcfpaths.get(0), e.getMessage());
                    System.exit(1);
                }
                ++n;
            }
        }
        return result;
    }
}
