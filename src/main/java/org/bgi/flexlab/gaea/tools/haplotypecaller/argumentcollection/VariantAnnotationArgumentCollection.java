package org.bgi.flexlab.gaea.tools.haplotypecaller.argumentcollection;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.bgi.flexlab.gaea.util.Utils;

public class VariantAnnotationArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Argument collection constructor which defines tool default annotation arguments when they are not overridden by the user
     *
     * @param defaultGroups                      List annotation group names to be used by default
     * @param defaultAnnotations                 List of annotation class names to be used by default
     * @param defaultAnnotationsToExclude        List of annotation class names to exclude by default. These override the default annotations and annotation groups.
     */
    public VariantAnnotationArgumentCollection(List<String> defaultGroups, List<String> defaultAnnotations, List<String> defaultAnnotationsToExclude) {
        Utils.nonNull(defaultGroups);
        Utils.nonNull(defaultAnnotations);
        Utils.nonNull(defaultAnnotationsToExclude);

        annotationGroupsToUse = new ArrayList<>(defaultGroups);
        annotationsToUse = new ArrayList<>(defaultAnnotations);
        annotationsToExclude = new ArrayList<>(defaultAnnotationsToExclude);
    }

    /**
     * Which annotations to include in variant calls in the output. These supplement annotations provided by annotation groups.
     */
    public List<String> annotationsToUse = new ArrayList<>();

    /**
     * Which annotations to exclude from output in the variant calls.  Note that this argument has higher priority than the
     * -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other
     * options.
     */
    public List<String> annotationsToExclude = new ArrayList<>();

    /**
     * Which groups of annotations to add to the output variant calls.
     * Any requirements that are not met (e.g. failing to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    public List<String> annotationGroupsToUse = new ArrayList<>();

}
