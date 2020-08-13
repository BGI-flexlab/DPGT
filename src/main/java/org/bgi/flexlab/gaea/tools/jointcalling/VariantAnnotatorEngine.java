package org.bgi.flexlab.gaea.tools.jointcalling;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.exception.UserException.BadArgumentValueException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.tools.haplotypecaller.annotator.VariantOverlapAnnotator;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.ReducibleAnnotationData;
import org.bgi.flexlab.gaea.tools.haplotypecaller.utils.RefMetaDataTracker;
import org.bgi.flexlab.gaea.tools.jointcalling.annotator.AnnotationInterfaceManager;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.GenotypeAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.InfoFieldAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.ReducibleAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.genotypegvcfs.annotation.VariantAnnotatorAnnotation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.GaeaGvcfVariantContextUtils;
import org.bgi.flexlab.gaea.tools.jointcalling.util.RodBinding;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VariantAnnotatorEngine {
	private List<InfoFieldAnnotation> requestedInfoAnnotations = Collections.emptyList();
	private List<InfoFieldAnnotation> requestedReducibleInfoAnnotations = new ArrayList<>();
	private List<InfoFieldAnnotation> requestedNonReducibleInfoAnnotations = new ArrayList<>();
	private List<GenotypeAnnotation> requestedGenotypeAnnotations = Collections.emptyList();
	private List<VAExpression> requestedExpressions = new ArrayList<>();
	private boolean expressionAlleleConcordance = false;
	private Set<String> overlapNames = new HashSet<String>();

	VariantOverlapAnnotator variantOverlapAnnotator = null;

	// Map of info field name to info field
	private final Map<String, VCFInfoHeaderLine> hInfoMap = new HashMap<>();

	protected static class VAExpression {

		public String fullName, fieldName;
		public RodBinding<VariantContext> binding;

		public VAExpression(String fullExpression, List<RodBinding<VariantContext>> bindings) {
			final int indexOfDot = fullExpression.lastIndexOf(".");
			if (indexOfDot == -1)
				throw new BadArgumentValueException(fullExpression, "it should be in rodname.value format");

			fullName = fullExpression;
			fieldName = fullExpression.substring(indexOfDot + 1);

			final String bindingName = fullExpression.substring(0, indexOfDot);
			for (final RodBinding<VariantContext> rod : bindings) {
				if (rod.getName().equals(bindingName)) {
					binding = rod;
					break;
				}
			}
		}
	}

	// use this constructor if you want to select specific annotations (and/or
	// interfaces)
	public VariantAnnotatorEngine(List<String> annotationGroupsToUse, List<String> annotationsToUse,
			List<String> annotationsToExclude) {
		initializeAnnotations(annotationGroupsToUse, annotationsToUse, annotationsToExclude);
		setReducibleAnnotations();
	}

	private void initializeAnnotations(List<String> annotationGroupsToUse, List<String> annotationsToUse,
			List<String> annotationsToExclude) {
		AnnotationInterfaceManager.validateAnnotations(annotationGroupsToUse, annotationsToUse);
		requestedInfoAnnotations = AnnotationInterfaceManager.createInfoFieldAnnotations(annotationGroupsToUse,
				annotationsToUse);
		requestedGenotypeAnnotations = AnnotationInterfaceManager.createGenotypeAnnotations(annotationGroupsToUse,
				annotationsToUse);
		excludeAnnotations(annotationsToExclude);
	}

	private void excludeAnnotations(List<String> annotationsToExclude) {
		if (annotationsToExclude.isEmpty())
			return;

		final List<InfoFieldAnnotation> tempRequestedInfoAnnotations = new ArrayList<>(requestedInfoAnnotations.size());
		for (final InfoFieldAnnotation annotation : requestedInfoAnnotations) {
			if (!annotationsToExclude.contains(annotation.getClass().getSimpleName()))
				tempRequestedInfoAnnotations.add(annotation);
		}
		requestedInfoAnnotations = tempRequestedInfoAnnotations;

		final List<GenotypeAnnotation> tempRequestedGenotypeAnnotations = new ArrayList<>(
				requestedGenotypeAnnotations.size());
		for (final GenotypeAnnotation annotation : requestedGenotypeAnnotations) {
			if (!annotationsToExclude.contains(annotation.getClass().getSimpleName()))
				tempRequestedGenotypeAnnotations.add(annotation);
		}
		requestedGenotypeAnnotations = tempRequestedGenotypeAnnotations;
	}
	
	public void initializeDBs(boolean hasDBSNP){
		if(hasDBSNP)
			overlapNames.add(VCFConstants.DBSNP_KEY);
	}

	public void initializeDBs(final ArrayList<VariantContext> dbSNPs,GenomeLocationParser parser) {
		// check to see whether comp rods were included
		ArrayList<VariantContext> dbSNPBindings = dbSNPs;
		if (dbSNPBindings != null && dbSNPBindings.size() == 0)
			dbSNPBindings = null;

		final Map<String, List<VariantContext>> overlapBindings = new LinkedHashMap<>();

		if (dbSNPBindings != null )
			overlapBindings.put( VCFConstants.DBSNP_KEY, dbSNPBindings); 

		variantOverlapAnnotator = new VariantOverlapAnnotator(dbSNPBindings, overlapBindings,parser);
	}

	public void invokeAnnotationInitializationMethods(final Set<VCFHeaderLine> headerLines,Set<String> sampleList) {
		for (final VariantAnnotatorAnnotation annotation : requestedInfoAnnotations) {
			annotation.initialize(headerLines,sampleList);
		}

		for (final VariantAnnotatorAnnotation annotation : requestedGenotypeAnnotations) {
			annotation.initialize(headerLines,sampleList);
		}
	}

	public Set<VCFHeaderLine> getVCFAnnotationDescriptions() {
		final Set<VCFHeaderLine> descriptions = new HashSet<>();

		for (final InfoFieldAnnotation annotation : requestedInfoAnnotations)
			descriptions.addAll(annotation.getDescriptions());
		for (final GenotypeAnnotation annotation : requestedGenotypeAnnotations)
			descriptions.addAll(annotation.getDescriptions());
		for (final String db : overlapNames) {
			if (VCFStandardHeaderLines.getInfoLine(db, false) != null)
				descriptions.add(VCFStandardHeaderLines.getInfoLine(db));
			else
				descriptions.add(new VCFInfoHeaderLine(db, 0, VCFHeaderLineType.Flag, db + " Membership"));
		}

		return descriptions;
	}

	private GenotypesContext annotateGenotypes(final RefMetaDataTracker tracker, final ChromosomeInformationShare ref, final VariantContext vc) {
		if (requestedGenotypeAnnotations.isEmpty())
			return vc.getGenotypes();

		final GenotypesContext genotypes = GenotypesContext.create(vc.getNSamples());
		for (final Genotype genotype : vc.getGenotypes()) {

			final GenotypeBuilder gb = new GenotypeBuilder(genotype);
			for (final GenotypeAnnotation annotation : requestedGenotypeAnnotations) {
				annotation.annotate(tracker, ref,  vc, genotype, gb);
			}
			genotypes.add(gb.make());
		}

		return genotypes;
	}

	public VariantContext annotateContext(final RefMetaDataTracker tracker, final ChromosomeInformationShare ref,
			final VariantContext vc) {
		// annotate genotypes
		final VariantContextBuilder builder = new VariantContextBuilder(vc)
				.genotypes(annotateGenotypes(tracker, ref, vc));
		VariantContext newGenotypeAnnotatedVC = builder.make();

		// annotate expressions where available
		final Map<String, Object> infoAnnotations = new LinkedHashMap<>(newGenotypeAnnotatedVC.getAttributes());
		annotateExpressions(tracker, tracker.getLocation(), newGenotypeAnnotatedVC, infoAnnotations);

		// go through all the requested info annotationTypes
		for (final InfoFieldAnnotation annotationType : requestedInfoAnnotations) {
			final Map<String, Object> annotationsFromCurrentType = annotationType.annotate(tracker, ref,
					newGenotypeAnnotatedVC);
			if (annotationsFromCurrentType != null)
				infoAnnotations.putAll(annotationsFromCurrentType);
		}

		// create a new VC in the with info and genotype annotations
		final VariantContext annotated = builder.attributes(infoAnnotations).make();

		// annotate db occurrences
		return annotateDBs( annotated);
	}
	
    private VariantContext annotateDBs(VariantContext vc) {
		Utils.nonNull(vc, "variantContext must not be null");
        return variantOverlapAnnotator.annotateOverlaps( variantOverlapAnnotator.annotateRsID(vc));
    }

	/**
	 * Combine (raw) data for reducible annotations (those that use raw data in
	 * gVCFs) Mutates annotationMap by removing the annotations that were
	 * combined
	 * 
	 * @param allelesList
	 *            the list of merged alleles across all variants being combined
	 * @param annotationMap
	 *            attributes of merged variant contexts -- is modifying by
	 *            removing successfully combined annotations
	 * @return a map containing the keys and raw values for the combined
	 *         annotations
	 */
	public Map<String, Object> combineAnnotations(final List<Allele> allelesList,
			Map<String, List<ReducibleAnnotationData>> annotationMap) {
		Map<String, Object> combinedAnnotations = new HashMap<>();

		// go through all the requested reducible info annotationTypes
		for (final InfoFieldAnnotation annotationType : requestedReducibleInfoAnnotations) {
			ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;
			if (annotationMap.containsKey(currentASannotation.getRawKeyName())) {
				final List<ReducibleAnnotationData> annotationValue = annotationMap
						.get(currentASannotation.getRawKeyName());
				final Map<String, Object> annotationsFromCurrentType = currentASannotation.combineRawData(allelesList,
						annotationValue);
				if(annotationsFromCurrentType == null)
					continue;
				combinedAnnotations.putAll(annotationsFromCurrentType);
				// remove the combined annotations so that the next method only
				// processes the non-reducible ones
				annotationMap.remove(currentASannotation.getRawKeyName());
			}
		}
		return combinedAnnotations;
	}

	/**
	 * Finalize reducible annotations (those that use raw data in gVCFs)
	 * 
	 * @param vc
	 *            the merged VC with the final set of alleles, possibly subset
	 *            to the number of maxAltAlleles for genotyping
	 * @param originalVC
	 *            the merged but non-subset VC that contains the full list of
	 *            merged alleles
	 * @return a VariantContext with the final annotation values for reducible
	 *         annotations
	 */
	public VariantContext finalizeAnnotations(VariantContext vc, VariantContext originalVC) {
		final Map<String, Object> infoAnnotations = new LinkedHashMap<>(vc.getAttributes());

		// go through all the requested info annotationTypes
		for (final InfoFieldAnnotation annotationType : requestedReducibleInfoAnnotations) {

			ReducibleAnnotation currentASannotation = (ReducibleAnnotation) annotationType;

			final Map<String, Object> annotationsFromCurrentType = currentASannotation.finalizeRawData(vc, originalVC);
			if (annotationsFromCurrentType != null) {
				infoAnnotations.putAll(annotationsFromCurrentType);
			}
			// clean up raw annotation data after annotations are finalized
			infoAnnotations.remove(currentASannotation.getRawKeyName());
		}

		// generate a new annotated VC
		final VariantContextBuilder builder = new VariantContextBuilder(vc).attributes(infoAnnotations);

		// annotate genotypes, creating another new VC in the process
		final VariantContext annotated = builder.make();
		return annotated;
	}

	private void setReducibleAnnotations() {
		for (final InfoFieldAnnotation annotationType : requestedInfoAnnotations) {
			if (annotationType instanceof ReducibleAnnotation)
				requestedReducibleInfoAnnotations.add(annotationType);
			else
				requestedNonReducibleInfoAnnotations.add(annotationType);
		}
	}
	
    private void annotateExpressions(final RefMetaDataTracker tracker, final GenomeLocation loc, final VariantContext vc, final Map<String, Object> infoAnnotations){

		Utils.nonNull(tracker, "refMetaDataTracker must not be null");
		Utils.nonNull(loc, "GenomeLocation must not be null");
		Utils.nonNull(vc, "VariantContext must not be null");
        // each requested expression
        for ( final VAExpression expression : requestedExpressions ) {

            // get the variant contexts for all the expressions at the location
            final Collection<VariantContext> expressionVCs = tracker.getValues(expression.binding, loc);
            if ( expressionVCs.isEmpty() )
                continue;

            // get the expression's variant context
            final VariantContext expressionVC = expressionVCs.iterator().next();

            // special-case the ID field
            if ( expression.fieldName.equals("ID") ) {
                if ( expressionVC.hasID() )
                    infoAnnotations.put(expression.fullName, expressionVC.getID());
            } else if (expression.fieldName.equals("ALT")) {
                infoAnnotations.put(expression.fullName, expressionVC.getAlternateAllele(0).getDisplayString());
            } else if (expression.fieldName.equals("FILTER")) {
                if ( expressionVC.isFiltered() ) {
                    infoAnnotations.put(expression.fullName, expressionVC.getFilters().toString().replace("[", "").replace("]", "").replace(" ", ""));
                } else {
                    infoAnnotations.put(expression.fullName, "PASS");
                }
            } else if ( expressionVC.hasAttribute(expression.fieldName) ) {
                // find the info field
                final VCFInfoHeaderLine hInfo = hInfoMap.get(expression.fullName);
                if ( hInfo == null ){
                    throw new UserException("Cannot annotate expression " + expression.fullName + " at " + loc + " for variant allele(s) " + vc.getAlleles() + ", missing header info");
                }

                //
                // Add the info field annotations
                //
                final boolean useRefAndAltAlleles = VCFHeaderLineCount.R == hInfo.getCountType();
                final boolean useAltAlleles = VCFHeaderLineCount.A == hInfo.getCountType();

                // Annotation uses ref and/or alt alleles or enforce allele concordance
                if ( (useAltAlleles || useRefAndAltAlleles) || expressionAlleleConcordance ){

                    // remove brackets and spaces from expression value
                    final String cleanedExpressionValue = expressionVC.getAttribute(expression.fieldName).toString().replaceAll("[\\[\\]\\s]", "");

                    // get comma separated expression values
                    final ArrayList<String> expressionValuesList = new ArrayList<String>(Arrays.asList(cleanedExpressionValue.split(",")));

                    // get the minimum biallelics without genotypes
                    final List<VariantContext> minBiallelicVCs = getMinRepresentationBiallelics(vc);
                    final List<VariantContext> minBiallelicExprVCs = getMinRepresentationBiallelics(expressionVC);

                    // check concordance
                    final List<String> annotationValues = new ArrayList<>();
                    boolean canAnnotate = false;
                    for ( final VariantContext biallelicVC : minBiallelicVCs ) {
                        // check that ref and alt alleles are the same
                        List<Allele> exprAlleles = biallelicVC.getAlleles();
                        boolean isAlleleConcordant = false;
                        int i = 0;
                        for ( final VariantContext biallelicExprVC : minBiallelicExprVCs ){
                            List<Allele> alleles = biallelicExprVC.getAlleles();
                            // concordant
                            if ( alleles.equals(exprAlleles) ){
                                // get the value for the reference if needed.
                                if ( i == 0 && useRefAndAltAlleles )
                                    annotationValues.add(expressionValuesList.get(i++));
                                // use annotation expression and add to vc
                                annotationValues.add(expressionValuesList.get(i));
                                isAlleleConcordant = true;
                                canAnnotate = true;
                                break;
                            }
                            i++;
                        }

                        // can not find allele match so set to annotation value to zero
                        if ( !isAlleleConcordant )
                            annotationValues.add("0");
                    }

                    // no allele matches so can not annotate
                    if ( !canAnnotate )
                        continue;

                    // add the annotation values
                    infoAnnotations.put(expression.fullName, annotationValues);
                } else {
                    // use all of the expression values
                    infoAnnotations.put(expression.fullName, expressionVC.getAttribute(expression.fieldName));
                }
            }
        }
    }
	
	private List<VariantContext> getMinRepresentationBiallelics(final VariantContext vc) {
        final List<VariantContext> minRepresentationBiallelicVCs = new ArrayList<VariantContext>();
        final boolean isMultiAllelic = vc.getNAlleles() > 2;
        if (isMultiAllelic) {
            final List<VariantContext> vcList = GaeaGvcfVariantContextUtils.splitVariantContextToBiallelics(vc);
            for (final VariantContext biallelicVC : vcList) {
                if (!biallelicVC.isSNP())
                    minRepresentationBiallelicVCs.add(GaeaGvcfVariantContextUtils.trimAlleles(biallelicVC, true, true));
                else
                    minRepresentationBiallelicVCs.add(biallelicVC);
            }
        } else {
            minRepresentationBiallelicVCs.add(vc);
        }

        return minRepresentationBiallelicVCs;
    }
	
	public List<InfoFieldAnnotation> getRequestedInfoAnnotations() {
        return requestedInfoAnnotations;
    }
}
