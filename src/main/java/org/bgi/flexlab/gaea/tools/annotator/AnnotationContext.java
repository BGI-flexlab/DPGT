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
 * Copyright (C)  2016  Pablo Cingolani(pcingola@users.sourceforge.net)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.ArrayUtils;
import org.bgi.flexlab.gaea.tools.annotator.effect.EffectType;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect;
import org.bgi.flexlab.gaea.tools.annotator.effect.VariantEffect.FunctionalClass;
import org.bgi.flexlab.gaea.tools.annotator.interval.*;

import java.io.Serializable;
import java.util.*;

/**
 * Common annotation context
 */
public class AnnotationContext implements Serializable{

	private static final long serialVersionUID = 5318310164237536996L;
	
	public static final String EFFECT_TYPE_SEPARATOR = "&";
	public static final String EFFECT_TYPE_SEPARATOR_SPLIT = "\\&";

	public static boolean debug = false;
	
	
	public static String ANN_FIELD_NAMES[] = { //
			"ALLELE", "GT", "GENOTYPE", //
			"EFFECT", "ANNOTATION", //
			"IMPACT", //
			"GENE", //
			"GENEID", //
			"FEATURE", //
			"FEATUREID", "TRID", //
			"BIOTYPE", //
			"RANK", "EXID", //
			"HGVS_C", "HGVS_DNA", "CODON", //
			"HGVS", "HGVS_P", "HGVS_PROT", "AA", //
			"HGVS_IVS", //
			"HGVS_Old", //
			"POS_CDNA", "CDNA_POS", //
			"LEN_CDNA", "CDNA_LEN", //
			"POS_CDS", "CDS_POS", //
			"LEN_CDS", "CDS_LEN", //
			"POS_AA", "AA_POS", //
			"LEN_AA", "AA_LEN", //
			"DISTANCE", //
			"VarType", //
			"Strand", //
			"ERRORS", "WARNINGS", "INFOS", //
	};


	Variant variant;
	Allele allele;
	String vcfFieldString; // Original 'raw' string from VCF Info field
	String vcfFieldStrings[]; // Original 'raw' strings from VCF info field: effectString.split()
	EffectType effectType;
	String effectTypesStr;
	List<EffectType> effectTypes;
	String effectDetails;
	int aaLen, aaPos;
	int cdsLen, cdsPos;
	int cDnaLen, cDnaPos;
	int distance;
	int rank, rankMax;
	BioType bioType;
	String codon, aa, hgvsC, hgvsP, hgvsProt, hgvsIVS, hgvsOld;
	VariantEffect.Coding coding;
	String genotype;
	String errorsWarnings;
	String geneName, geneId, featureType, featureId, transcriptId, exonId;
	VariantEffect.EffectImpact impact;
	VariantEffect.FunctionalClass funClass;
	VariantEffect variantEffect;
	boolean useSequenceOntology;
	boolean useSimpleEffectNom;
	boolean useHgvs;
	boolean useGeneId;
	boolean useFirstEffect;
	String strand;

	private Map<String, Object> annoItems;

	public AnnotationContext() {
		init();
	}

	public AnnotationContext(VariantEffect variantEffect) {
		this(variantEffect, true, false);
	}

	public AnnotationContext(VariantEffect variantEffect, boolean useSequenceOntology, boolean useFirstEffect) {
		init();
		this.variantEffect = variantEffect;
		this.useSequenceOntology = useSequenceOntology;
		this.useFirstEffect = useFirstEffect;
		set(variantEffect);
	}

	public AnnotationContext(VariantEffect variantEffect, boolean useSimpleEffectNom) {
		init();
		this.variantEffect = variantEffect;
		this.useSimpleEffectNom = useSimpleEffectNom;
		this.useSequenceOntology = true;
		this.useFirstEffect = false;
		set(variantEffect);
	}
	
	void init() {
		aaLen = aaPos = cdsLen = cdsPos = cDnaLen = cDnaPos = distance = rank = rankMax = -1;
		vcfFieldString = effectTypesStr = effectDetails = codon = aa = hgvsC = hgvsP = hgvsProt = genotype = errorsWarnings = geneName = geneId = featureType = featureId = transcriptId = exonId = errorsWarnings = "";
		hgvsIVS = hgvsOld = ".";
		bioType = null;
		impact = null;
		funClass = FunctionalClass.NONE;
		useSequenceOntology = true;
		useHgvs = true;
		useGeneId = false;
		strand = ".";
		annoItems = new HashMap<>();
	}


	 /**
     * @return the annoItem map
     */
    public Map<String, Object> getAnnoItems() {
        return Collections.unmodifiableMap(annoItems);
    }
    
    // todo -- define common annoItems as enum

    public void setAnnoItems(Map<String, ?> map) {
        clearAnnoItems();
        putAnnoItems(map);
    }

	public void putAnnoItems(Map<String, ?> map) {
        if ( map != null ) {
            if (annoItems.isEmpty()) {
                annoItems.putAll(map);
            } else {
                for ( Map.Entry<String, ?> elt : map.entrySet() ) {
                	putAnnoItem(elt.getKey(), elt.getValue(), false);
                }
            }
        }
    }
	
	public void putAnnoItem(String key, Object value) {
		putAnnoItem(key, value, false);
	}

	public void putAnnoItem(String key, Object value, boolean allowOverwrites) {
        if ( ! allowOverwrites && hasAnnoItem(key) )
            throw new IllegalStateException("Attempting to overwrite key->value binding: key = " + key + " this = " + this);
        annoItems.put(key, value);
    }

	
	 public boolean hasAnnoItem(String key) {
		 return annoItems.containsKey(key);
	}


	public void clearAnnoItems() {
		annoItems = new HashMap<>();
	}
	
    /**
     * @param key    the AnnoItem key
     *
     * @return the AnnoItem value for the given key (or null if not set)
     */
    public Object getAnnoItem(String key) {
    	if(hasAnnoItem(key))
			return annoItems.get(key);
		if(ArrayUtils.contains(ANN_FIELD_NAMES, key))
			return getFieldByName(key);
        return null;
    }

    public Object getAnnoItem(String key, Object defaultValue) {
        if ( hasAnnoItem(key) )
            return annoItems.get(key);
        else if(ArrayUtils.contains(ANN_FIELD_NAMES, key))
			return getFieldByName(key);
        else
            return defaultValue;
    }

    /** returns the value as an empty list if the key was not found,
        as a java.util.List if the value is a List or an Array,
        as a Collections.singletonList if there is only one value */
    @SuppressWarnings("unchecked")
    public List<Object> getAnnoItemAsList(String key) {
        Object o = getAnnoItem(key);
        if ( o == null ) return Collections.emptyList();
        if ( o instanceof List ) return (List<Object>)o;
        if ( o.getClass().isArray() ) return Arrays.asList((Object[])o);
        return Collections.singletonList(o);
    }

    public String getAnnoItemAsString(String key, String defaultValue) {
        Object x = getAnnoItem(key);
        if ( x == null || x.equals("")) return defaultValue;
        if ( x instanceof String ) return (String)x;
        return String.valueOf(x); // throws an exception if this isn't a string
    }

    public int getAnnoItemAsInt(String key, int defaultValue) {
        Object x = getAnnoItem(key);
        if ( x == null || x == VCFConstants.MISSING_VALUE_v4 ) return defaultValue;
        if ( x instanceof Integer ) return (Integer)x;
        return Integer.valueOf((String)x); // throws an exception if this isn't a string
    }

    public double getAnnoItemAsDouble(String key, double defaultValue) {
        Object x = getAnnoItem(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Double ) return (Double)x;
        if ( x instanceof Integer ) return (Integer)x;
        return Double.valueOf((String)x); // throws an exception if this isn't a string
    }

    public boolean getAnnoItemAsBoolean(String key, boolean defaultValue) {
        Object x = getAnnoItem(key);
        if ( x == null ) return defaultValue;
        if ( x instanceof Boolean ) return (Boolean)x;
        return Boolean.valueOf((String)x); // throws an exception if this isn't a string
    }

	public void addEffectType(EffectType effectType) {
		effectTypes.add(effectType);
		this.effectType = null;
	}


	public String getAa() {
		return aa;
	}

	public int getAaLen() {
		return aaLen;
	}

	public int getAaPos() {
		return aaPos;
	}

	public String getAlt() {
		return getGenotype();
	}

	public BioType getBioType() {
		return bioType;
	}

	public int getcDnaLen() {
		return cDnaLen;
	}

	public int getcDnaPos() {
		return cDnaPos;
	}

	public int getCdsLen() {
		return cdsLen;
	}

	public int getCdsPos() {
		return cdsPos;
	}

	public VariantEffect.Coding getCoding() {
		return coding;
	}

	public String getCodon() {
		return codon;
	}

	public int getDistance() {
		return distance;
	}

	public String getEffectDetails() {
		return effectDetails;
	}

	public String getEffectsStr() {
		StringBuilder sb = new StringBuilder();
		for (EffectType et : effectTypes) {
			if (sb.length() > 0) sb.append(EFFECT_TYPE_SEPARATOR);
			sb.append(et);
		}
		return sb.toString();
	}

	public String getEffectsStrSo() {
		StringBuilder sb = new StringBuilder();
		for (EffectType et : effectTypes) {
			if (sb.length() > 0) sb.append(EFFECT_TYPE_SEPARATOR);
			sb.append(et.toSequenceOntology(null));
		}
		return sb.toString();
	}

	public EffectType getEffectType() {
		if (effectType != null) return effectType;
		if (effectTypes == null || effectTypes.isEmpty()) return EffectType.NONE;

		// Pick highest effect type
		effectType = EffectType.NONE;
		for (EffectType et : effectTypes)
			if (et.compareTo(effectType) < 0) effectType = et;

		return effectType;
	}

	public List<EffectType> getEffectTypes() {
		return effectTypes;
	}

	public String getEffectTypesStr() {
		return effectTypesStr;
	}

	public String getErrorsWarning() {
		return errorsWarnings;
	}

	public String getExonId() {
		return exonId;
	}

	public String getFeatureId() {
		return featureId;
	}

	public String getFeatureType() {
		return featureType;
	}
	
	/**
	 * Get a subfield by name
	 */
	public String getFieldByName(String fieldName) {
		switch (fieldName) {

		case "ALLELE":
		case "GT":
		case "GENOTYPE":
		case "GENOTYPE_NUMBER":
			return genotype;
			

		case "EFFECT":
		case "ANNOTATION":
			return effectTypesStr;

		case "IMPACT":
			return impact != null ? impact.toString() : ".";

		case "FUNCLASS":
			return funClass != null ? funClass.toString() : ".";

		case "GENE":
			return geneName;

		case "GENEID":
			return geneId;

		case "FEATURE":
		case "FEATURE_TYPE":
			return featureType;

		case "FEATUREID":
			return featureId;

		case "TRID":
			return transcriptId;

		case "BIOTYPE":
			return (bioType == null ? "." : bioType.toString());

		case "RANK":
		case "EXON_RANK":
			return Integer.toString(rank);

		case "EXID":
		case "EXON_ID":
			return exonId;

		case "RANK_MAX":
			return Integer.toString(rankMax);

		case "HGVS_C":
		case "HGVS_DNA":
			return hgvsC;

		case "HGVS_IVS":
			return hgvsIVS;

		case "HGVS_Old":
			return hgvsOld;

		case "CODON":
			return codon;

		case "HGVS":
		case "HGVS_P":
			return hgvsP;

		case "HGVS_PROT":
			return hgvsProt;

		case "AA":
			return aa;

		case "POS_CDNA":
		case "CDNA_POS":
			return Integer.toString(cDnaPos);

		case "LEN_CDNA":
		case "CDNA_LEN":
			return Integer.toString(cDnaLen);

		case "POS_CDS":
		case "CDS_POS":
			return Integer.toString(cdsPos);

		case "LEN_CDS":
		case "CDS_LEN":
			return Integer.toString(cdsLen);

		case "POS_AA":
		case "AA_POS":
			return Integer.toString(aaPos);

		case "LEN_AA":
		case "AA_LEN":
			return Integer.toString(aaLen);

		case "CODING":
			return coding != null ? coding.toString() : ".";

		case "DISTANCE":
			return Integer.toString(distance);

		case "VarType":
			return variant.getVariantType().toString();

		case "Strand":
			return strand;

		case "ERRORS":
		case "WARNINGS":
		case "INFOS":
			return errorsWarnings;

		default:
			throw new RuntimeException("Field '" + fieldName + "' not found.");
		}
	}


	public VariantEffect.FunctionalClass getFunClass() {
		return funClass;
	}

	public String getGeneId() {
		return geneId;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getGenotype() {
		return genotype;
	}

	public String getHgvsC() {
		return hgvsC;
	}

	public String getHgvsDna() {
		return hgvsC;
	}

	public String getHgvsP() {
		return hgvsP;
	}

	public String getHgvsProt() {
		return hgvsP;
	}

	public VariantEffect.EffectImpact getImpact() {
		return impact;
	}

	public int getRank() {
		return rank;
	}

	public int getRankMax() {
		return rankMax;
	}

	public String getTranscriptId() {
		return transcriptId;
	}

	/**
	 * String from VCF file (original, unparsed, string)
	 */
	public String getVcfFieldString() {
		return vcfFieldString;
	}

	/**
	 * Get a subfield as an index
	 */
	public String getVcfFieldString(int index) {
		if (index >= vcfFieldStrings.length) return null;
		return vcfFieldStrings[index];
	}

	
	/**
	 * Set all fields from 'variantEffect'
	 */
	void set(VariantEffect variantEffect) {
		// Allele
		variant = variantEffect.getVariant();
		Gene gene = variantEffect.getGene();
		Marker marker = variantEffect.getMarker();
		Transcript tr = variantEffect.getTranscript();

		// Genotype
		if (variant.getGenotype() != null) genotype = variant.getGenotype();
		else if (!variant.isVariant()) genotype = variant.getReference();
		else genotype = variant.getAlt();
//		 else if (var.isNonRef()) genotype = var.getGenotype();

//		System.err.println("---2"+variant.getAlt());
//		//TODO  for del: CT->C
		// 		variant.getAlt() : del: T->
//		//      variant.getGenotype() : del: CT->C
		allele = Allele.create(variant.getGenotype(), false);

		// Effect
		effectType = variantEffect.getEffectType();
		effectTypes = variantEffect.getEffectTypes();

		effectTypesStr = variantEffect.getEffectTypeString(true, useFirstEffect);
		effectTypesStr = variantEffect.getSimpleEffectTypeString(useFirstEffect);

		// Impact
		impact = variantEffect.getEffectImpact();

		// Functional class
		funClass = variantEffect.getFunctionalClass();

		// Gene
		if (gene != null) {
			if (variantEffect.isMultipleGenes()) {
				setGeneNameIdMultiple(variantEffect);
			} else {
				geneName = gene.getGeneName();
				geneId = gene.getId();
			}
		} else if (marker instanceof Intergenic) {
			geneName = ((Intergenic) marker).getName();
			geneId = marker.getId();
		} else {
			geneName = geneId = "";
		}

		// Feature type & ID
		featureType = featureId = "";
		if (marker != null) {
			if (marker instanceof Custom) {
				// Custom
				featureType = marker.getType() + EFFECT_TYPE_SEPARATOR + ((Custom) marker).getLabel();
				featureId = marker.getId();
			} else if (marker instanceof Regulation) {
				// Regulation includes cell type
				Regulation reg = (Regulation) marker;
				featureType = reg.getType() + EFFECT_TYPE_SEPARATOR + reg.getName() + ":" + reg.getCellType();
				featureId = marker.getId();
			} else if (marker instanceof ProteinStructuralInteractionLocus) {
				featureType = "interaction";
				featureId = marker.getId();
			} else if (tr != null) {
				featureType = "transcript";
				featureId = transcriptId = tr.getId();
				strand = tr.getStrand();
				// Append version number (this is recommended by HGVS specification)
				if (tr.getVersion() != null && !tr.getVersion().isEmpty()) featureId += "." + tr.getVersion();
			} else {
				featureType = marker.getType().toSequenceOntology(null);
				featureId = marker.getId();
			}
		}

		// Biotype
		if (tr != null) {
			if ((tr.getBioType() != null) && (tr.getBioType() != null)) {
				bioType = tr.getBioType();
			} else {
				// No biotype? Add protein_coding of we know it is.
				bioType = BioType.coding(tr.isProteinCoding());
			}
		} else {
			bioType = null;
		}

		// Rank
		Exon ex = variantEffect.getExon();
		rank = -1;
		if (ex != null) {
			rank = ex.getRank();
			rankMax = tr.numChilds();
		} else {
			// Do we have an intron?
			Intron intron = variantEffect.getIntron();
			if (intron != null) {
				rank = intron.getRank();
				rankMax = Math.max(0, tr.numChilds() - 1);
			} else if (tr != null && marker != null) {
				// Can we try to find an exon?
				for (Exon e : tr)
					if (e.intersects(marker)) {
						rank = e.getRank();
						rankMax = tr.numChilds();
						break;
					}
			}
		}

		// Codon change
		codon = variantEffect.getCodonChangeMax();

		// AA change
		aa = variantEffect.getAaChange();

		// HGVS notation
		hgvsC = variantEffect.getHgvsDna();
		hgvsP = variantEffect.getHgvsP();
		hgvsProt = variantEffect.getHgvsProt();
		hgvsIVS = variantEffect.getHgvsIVS();
		hgvsOld = variantEffect.getHgvsOld();

		// cDna position & len (cDNA is the DNA version of mRNA)
		if (tr != null) {
			cDnaPos = variantEffect.getcDnaPos();
			if (cDnaPos >= 0 ) cDnaPos++; // 1-based position;
			cDnaLen = tr.mRna().length();
		} else {
			cDnaPos = cDnaLen = -1;
		}

		// CDS position / length
		if (tr != null) {
			cdsPos = variantEffect.getCodonNum() * 3 + variantEffect.getCodonIndex();
			if (cdsPos >= 0) cdsPos++; // 1-based position;
			cdsLen = variantEffect.getCdsLength();
		} else {
			cdsPos = cdsLen = -1;
		}

		// Protein position / protein length
		if (tr != null) {
			aaPos = variantEffect.getCodonNum();
			if (aaPos >= 0) aaPos++; // 1-based position;
			aaLen = variantEffect.getAaLength();
		} else {
			aaPos = aaLen = -1;
		}

		// Distance: Mostly used for non-coding variants
		distance = variantEffect.getDistance();

		if (variantEffect.hasError() || variantEffect.hasWarning()) {
			StringBuilder err = new StringBuilder();
			// Add errors
			if (!variantEffect.getError().isEmpty()) {
				err.append(variantEffect.getError());
			}

			// Add warnings
			if (!variantEffect.getWarning().isEmpty()) {
				if (err.length() > 0) err.append(EFFECT_TYPE_SEPARATOR);
				err.append(variantEffect.getWarning());
			}

			errorsWarnings = err.toString();
		}

	}

	public void setAa(String aa) {
		this.aa = aa;
	}

	public void setAaLen(int aaLen) {
		this.aaLen = aaLen;
	}

	public void setBioType(BioType bioType) {
		this.bioType = bioType;
	}

	public void setCoding(VariantEffect.Coding coding) {
		this.coding = coding;
	}

	public void setCodon(String codon) {
		this.codon = codon;
	}

	public void setEffectDetails(String effectDetails) {
		this.effectDetails = effectDetails;
	}

	public void setEffectType(EffectType effect) {
		effectTypes = new LinkedList<EffectType>();
		addEffectType(effect);
	}

	public void setExonId(String exonId) {
		this.exonId = exonId;
	}


	public void setFunClass(VariantEffect.FunctionalClass funClass) {
		this.funClass = funClass;
	}

	public void setGeneId(String geneId) {
		this.geneId = geneId;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	/**
	 * Structural variant having multiple genes: Set all geneNames and geneIds
	 */
	void setGeneNameIdMultiple(VariantEffect variantEffect) {
		// Get all genes
		List<Gene> genes = variantEffect.getGenes();

		// Sort by geneName
		Collections.sort(genes, new Comparator<Gene>() {
			@Override
			public int compare(Gene g1, Gene g2) {
				return g1.getGeneName().compareTo(g2.getGeneName());
			}
		});

		// Create gene name and geneId strings
		StringBuilder geneNames = new StringBuilder();
		StringBuilder geneIds = new StringBuilder();

		String sep = "";
		for (Gene g : genes) {
			geneNames.append(sep + g.getGeneName());
			geneIds.append(sep + g.getId());
			if (sep.isEmpty()) sep = EFFECT_TYPE_SEPARATOR;
		}

		geneName = geneNames.toString();
		geneId = geneIds.toString();
	}

	public void setGenotype(String genotype) {
		this.genotype = genotype;
	}

	public void setImpact(VariantEffect.EffectImpact impact) {
		this.impact = impact;
	}

	public void setTranscriptId(String transcriptId) {
		this.transcriptId = transcriptId;
	}

	public void setUseFirstEffect(boolean useFirstEffect) {
		this.useFirstEffect = useFirstEffect;
	}

	public void setUseGeneId(boolean useGeneId) {
		this.useGeneId = useGeneId;
	}

	public void setUseHgvs(boolean useHgvs) {
		this.useHgvs = useHgvs;
	}

	public Variant getVariant() {
		return variant;
	}

	public void setVariant(Variant variant) {
		this.variant = variant;
	}

	public void setAllele(Allele allele) {
		this.allele = allele;
	}

	public void setAllele(String alt) {
		this.allele = Allele.create(alt, false);
	}

	public Allele getAllele() {
		return allele;
	}

	public List<String> getAnnoValues(List<String> fields){
		List<String> annoValues = new ArrayList<>();
		for (String field : fields){
			annoValues.add(getAnnoItemAsString(field, "."));
		}
		return annoValues;
	}

	public String toAnnoString(List<String> fields){
		getAnnoValues(fields);
		return String.join("|", getAnnoValues(fields));
	}

}
