package org.bgi.flexlab.dpgt;

import java.nio.file.Paths;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Arrays;
import java.util.Set;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKVariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationPluginDescriptor;
import org.broadinstitute.barclay.argparser.ClassFinder;
import org.broadinstitute.barclay.argparser.CommandLineException;


public class GVCFsSyncGenotyper {
    private ReferenceDataSource reference;
    private MultiVariantSyncReader reader;
    VariantAnnotatorEngine annotationEngine;

    public GVCFsSyncGenotyper(final String refpath, final List<String> vcfpaths, final SimpleInterval interval) {
        // init reference
        reference = ReferenceDataSource.of(Paths.get(refpath));
        reader = new MultiVariantSyncReader();
        reader.open(vcfpaths);
        reader.query(interval);
        List<Annotation> annotations = makeVariantAnnotations();
        annotationEngine = new VariantAnnotatorEngine(annotations, null, Collections.emptyList(), false);
    }

    public void run() {
        ArrayList<VariantContext> variantContexts;
        while(!(variantContexts = reader.read()).isEmpty()) {
            final SimpleInterval variantInterval = new SimpleInterval(variantContexts.get(0));
            // apply(variantInterval, variantContexts, new ReferenceContext(reference, variantInterval)); 
        }
    }

    // public void apply(final Locatable loc, List<VariantContext> variants, ReferenceContext ref) {

    // }


    private List<Annotation> makeVariantAnnotations() {
        GATKAnnotationArgumentCollection userArgs = new DefaultGATKVariantAnnotationArgumentCollection();
        GATKAnnotationPluginDescriptor pluginDescriptor = new GATKAnnotationPluginDescriptor(userArgs, Collections.emptyList(), Arrays.asList(StandardAnnotation.class));
        findPluginsForDescriptor(pluginDescriptor);
        pluginDescriptor.validateAndResolvePlugins();
        return pluginDescriptor.getResolvedInstances();
    }

    private void findPluginsForDescriptor(
        GATKAnnotationPluginDescriptor pluginDescriptor) {
        final ClassFinder classFinder = new ClassFinder();
        pluginDescriptor.getPackageNames().forEach(
                pkg -> classFinder.find(pkg, pluginDescriptor.getPluginBaseClass()));
        final Set<Class<?>> pluginClasses = classFinder.getClasses();

        final List<Object> plugins = new ArrayList<>(pluginClasses.size());
        for (Class<?> c : pluginClasses) {
            if (pluginDescriptor.includePluginClass(c)) {
                try {
                    final Object plugin = pluginDescriptor.createInstanceForPlugin(c);
                    plugins.add(plugin);
                } catch (InstantiationException | IllegalAccessException e) {
                    throw new CommandLineException.CommandLineParserInternalException("Problem making an instance of plugin " + c +
                            " Do check that the class has a non-arg constructor", e);
                }
            }
        }
    }
}
