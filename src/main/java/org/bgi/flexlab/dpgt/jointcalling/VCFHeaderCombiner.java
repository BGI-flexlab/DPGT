package org.bgi.flexlab.dpgt.jointcalling;
import java.lang.String;

public class VCFHeaderCombiner {
    public native void Combine(String[] vcfpaths, String outpath);
}
