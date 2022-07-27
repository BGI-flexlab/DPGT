package org.bgi.flexlab.dpgt.jointcalling;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import java.util.BitSet;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class VariantSiteSet {
    private static final Logger logger = LoggerFactory.getLogger(VariantSiteSet.class);
    public SimpleInterval interval;
    public BitSet data;

    public VariantSiteSet(final SimpleInterval interval) {
        this.interval = interval;
        this.data = new BitSet(this.interval.size());
    }

    

}
