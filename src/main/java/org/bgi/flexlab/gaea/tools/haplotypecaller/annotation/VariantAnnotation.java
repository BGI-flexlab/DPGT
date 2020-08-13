package org.bgi.flexlab.gaea.tools.haplotypecaller.annotation;

import java.util.List;

public abstract class VariantAnnotation {

    /**
     * Return the keys
     */
    public abstract List<String> getKeyNames();

    @Override
    public String toString() {
        return getClass().getSimpleName();
    }
}