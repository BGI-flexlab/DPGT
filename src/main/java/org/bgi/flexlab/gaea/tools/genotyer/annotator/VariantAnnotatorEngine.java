/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

package org.bgi.flexlab.gaea.tools.genotyer.annotator;


import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.bgi.flexlab.gaea.data.structure.pileup.Mpileup;
import org.bgi.flexlab.gaea.data.structure.pileup.Pileup;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.vcf.VariantDataTracker;
import org.bgi.flexlab.gaea.tools.genotyer.annotator.interfaces.*;
import org.bgi.flexlab.gaea.tools.genotyer.genotypeLikelihoodCalculator.PerReadAlleleLikelihoodMap;

import java.util.*;

/**
 *  rm database part cause this part take lots of time for large scale dbsnp etc. We will annotate this at
 *  GaeaAnnotator which is a database version and much faster.
 */
public class VariantAnnotatorEngine {

    private List<InfoFieldAnnotation> requestedInfoAnnotations = Collections.emptyList();
    private List<GenotypeAnnotation> requestedGenotypeAnnotations = Collections.emptyList();
    //private List<VAExpression> requestedExpressions = new ArrayList<VAExpression>();

    private final HashMap<VariantDataTracker, String> dbAnnotations = new HashMap<VariantDataTracker, String>();
    //private final AnnotatorCompatible walker;

    private boolean requireStrictAlleleMatch = false;

    // use this constructor if you want all possible annotations
    public VariantAnnotatorEngine(List<String> annotationsToExclude) {
        requestedInfoAnnotations = AnnotationInterfaceManager.createAllInfoFieldAnnotations();
        requestedGenotypeAnnotations = AnnotationInterfaceManager.createAllGenotypeAnnotations();
        excludeAnnotations(annotationsToExclude);
    }

    // use this constructor if you want to select specific annotations (and/or interfaces)
    public VariantAnnotatorEngine(List<String> annotationGroupsToUse, List<String> annotationsToUse, List<String> annotationsToExclude) {
        initializeAnnotations(annotationGroupsToUse, annotationsToUse, annotationsToExclude);
    }

    // experimental constructor for active region traversal
    public VariantAnnotatorEngine() {
        requestedInfoAnnotations = AnnotationInterfaceManager.createInfoFieldAnnotations(Arrays.asList("ActiveRegionBasedAnnotation"), Collections.<String>emptyList());
    }

    private void initializeAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse, List<String> annotationsToExclude) {
        AnnotationInterfaceManager.validateAnnotations(annotationGroupsToUse, annotationsToUse);
        requestedInfoAnnotations = AnnotationInterfaceManager.createInfoFieldAnnotations(annotationGroupsToUse, annotationsToUse);
        requestedGenotypeAnnotations = AnnotationInterfaceManager.createGenotypeAnnotations(annotationGroupsToUse, annotationsToUse);
        excludeAnnotations(annotationsToExclude);
    }

    private void excludeAnnotations(List<String> annotationsToExclude) {
        if (annotationsToExclude == null || annotationsToExclude.size() == 0 )
            return;

        List<InfoFieldAnnotation> tempRequestedInfoAnnotations = new ArrayList<InfoFieldAnnotation>(requestedInfoAnnotations.size());
        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations ) {
            if ( !annotationsToExclude.contains(annotation.getClass().getSimpleName()) )
                tempRequestedInfoAnnotations.add(annotation);
        }
        requestedInfoAnnotations = tempRequestedInfoAnnotations;

        List<GenotypeAnnotation> tempRequestedGenotypeAnnotations = new ArrayList<GenotypeAnnotation>(requestedGenotypeAnnotations.size());
        for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
            if ( !annotationsToExclude.contains(annotation.getClass().getSimpleName()) )
                tempRequestedGenotypeAnnotations.add(annotation);
        }
        requestedGenotypeAnnotations = tempRequestedGenotypeAnnotations;
    }

    public void invokeAnnotationInitializationMethods( Set<VCFHeaderLine> headerLines ) {
        for ( VariantAnnotatorAnnotation annotation : requestedInfoAnnotations ) {
            annotation.initialize(headerLines);
        }

        for ( VariantAnnotatorAnnotation annotation : requestedGenotypeAnnotations ) {
            annotation.initialize(headerLines);
        }
    }

    public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {

        Set<VCFHeaderLine> descriptions = new HashSet<VCFHeaderLine>();
       // System.out.println("requestedInfoAnnotations :====");
       // System.out.println(requestedInfoAnnotations.size());
        for ( InfoFieldAnnotation annotation : requestedInfoAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( GenotypeAnnotation annotation : requestedGenotypeAnnotations )
            descriptions.addAll(annotation.getDescriptions());
        for ( String db : dbAnnotations.values() ) {
            if ( VCFStandardHeaderLines.getInfoLine(db, false) != null )
                descriptions.add(VCFStandardHeaderLines.getInfoLine(db));
            else
                descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, db + " Membership"));
        }

        return descriptions;
    }

    public void setRequireStrictAlleleMatch( final boolean requireStrictAlleleMatch ) {
        this.requireStrictAlleleMatch = requireStrictAlleleMatch;
    }

    public VariantContext annotateContext(
    									   final VariantDataTracker tracker,
                                           final ChromosomeInformationShare ref,
                                           final Mpileup mpileup,
                                           VariantContext vc) {
        return annotateContext(tracker,ref, mpileup, vc, null);
    }

    public VariantContext annotateContext(
    									  final VariantDataTracker tracker,
                                          final ChromosomeInformationShare ref,
                                          final Mpileup mpileup,
                                          VariantContext vc,
                                          final Map<String,PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap) {
        Map<String, Object> infoAnnotations = new LinkedHashMap<String, Object>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, ref, mpileup, vc, perReadAlleleLikelihoodMap);
            if ( annotationsFromCurrentType != null )
                infoAnnotations.putAll(annotationsFromCurrentType);
            /*else {
                System.err.println("variant:" + vc.getStart() + "-" + vc.getEnd() + " error annotation:" + annotationType.getKeyNames());
            }*/
        }
       // if(vc.getStart()==145413)
       // 	System.out.println("vc annotator3:"+vc.toString());
        
        // generate a new annotated VC
        VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);
       // if(vc.getStart()==145413)
       //	System.out.println("vc annotator4:"+vc.toString());
        
        // annotate genotypes, creating another new VC in the process
        return builder.genotypes(annotateGenotypes(tracker, ref, mpileup, vc, perReadAlleleLikelihoodMap)).make();
    }

    public VariantContext annotateContext(final Map<String, PerReadAlleleLikelihoodMap> perReadAlleleLikelihoodMap, VariantContext vc) {
        Map<String, Object> infoAnnotations = new LinkedHashMap<String, Object>(vc.getAttributes());

        // go through all the requested info annotationTypes
        for ( InfoFieldAnnotation annotationType : requestedInfoAnnotations ) {
            if ( !(annotationType instanceof ActiveRegionBasedAnnotation) )
                continue;

            Map<String, Object> annotationsFromCurrentType = annotationType.annotate(perReadAlleleLikelihoodMap, vc);
            if ( annotationsFromCurrentType != null ) {
                infoAnnotations.putAll(annotationsFromCurrentType);
            }
        }

        // generate a new annotated VC
        VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);

        // annotate genotypes, creating another new VC in the process
        return builder.genotypes(annotateGenotypes(null, null, null, vc, perReadAlleleLikelihoodMap)).make();
    }

    private GenotypesContext annotateGenotypes(final VariantDataTracker tracker,
                                               final ChromosomeInformationShare ref, final Mpileup mpileup,
                                               final VariantContext vc,
                                               final Map<String,PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( requestedGenotypeAnnotations.isEmpty() )
            return vc.getGenotypes();

        final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
        for ( final Genotype genotype : vc.getGenotypes() ) {
           // AlignmentContext context = null;
            PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap = null;
          //  if (mpileup != null)
          //      context = stratifiedContexts.get(genotype.getSampleName());
            if (stratifiedPerReadAlleleLikelihoodMap != null)
                perReadAlleleLikelihoodMap = stratifiedPerReadAlleleLikelihoodMap.get(genotype.getSampleName());

           // System.out.println("sample:" + genotype.getSampleName());
            Pileup pileup = mpileup.getCurrentPosPileup().get(genotype.getSampleName());
            final GenotypeBuilder gb = new GenotypeBuilder(genotype);
            for ( final GenotypeAnnotation annotation : requestedGenotypeAnnotations ) {
                annotation.annotate(tracker, ref, pileup, vc, genotype, gb, perReadAlleleLikelihoodMap);
            }
            genotypes.add(gb.make());
        }

        return genotypes;
    }
}
