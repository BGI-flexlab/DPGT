package org.bgi.flexlab.gaea.data.structure.variant.statistic;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;

import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFCodec;

public class VariantBasicStatistic {
	// static array
	private IntArray basicArray = null;

	private static final int BASIC_LENGTH = VariantEnum.TRANSVERSIONS.ordinal() + 1;

	private static final int ALL_LENGTH = VariantEnum.values().length;

	private static final String mode1 = "%-40s\t%-40S\n";

	private String sampleName = null;

	private HashMap<String, IntArray> sampleStatics = new HashMap<String, IntArray>();

	private boolean individuals = true;
	
	private int mixed = 0;

	public VariantBasicStatistic(boolean individuals) {
		this.individuals = individuals;
		if (individuals)
			basicArray = new IntArray(ALL_LENGTH, 0);
		else
			basicArray = new IntArray(BASIC_LENGTH, 0);

		sampleStatics.clear();
	}

	public void variantStatic(VariantContext vc) {
		Allele ref = vc.getReference();
		List<Allele> alt = vc.getAlternateAlleles();
		boolean snp = false, mnp = false, insertion = false, deletion = false;

		if(vc.isVariant())
			basicArray.incr(VariantEnum.VARIANTS);
		if (vc.isSNP())
			snp = true;
		else if (vc.isMNP())
			mnp = true;
		else if (vc.isIndel()) {
			for (Allele allele : alt) {
				if (ref.length() > allele.length())
					deletion = true;
				if (ref.length() < allele.length())
					insertion = true;
			}
		} else if (vc.isMixed()) {
			for (Allele allele : alt) {
				if (ref.length() == allele.length()) {
					if (ref.length() == 1)
						snp = true;
					else
						mnp = true;
				} else {
					if (ref.length() > allele.length())
						deletion = true;
					if (ref.length() < allele.length())
						insertion = true;
				}
			}
		}

		boolean[] conf = new boolean[] { snp, mnp, insertion, deletion };

		basicStatics(vc, basicArray, conf);

		if (basicArray.length() != BASIC_LENGTH)
			genotyStatic(vc, conf);
	}

	private void basicStatics(VariantContext vc, IntArray array, boolean[] conf) {
		if (conf[0])
			array.incr(VariantEnum.SNPS);
		if (conf[1])
			array.incr(VariantEnum.MNPS);
		if (conf[2]) {
			array.incr(VariantEnum.INSERTIONS);
		}
		if (conf[3]) {
			array.incr(VariantEnum.DELETIONS);
		}

		if (conf[2] || conf[3])
			array.incr(VariantEnum.INDELS);

		if (conf[0] && vc.isBiallelic() && VariantContextUtils.isTransition(vc)) {
			array.incr(VariantEnum.TRANSITIONS);
		} else if (conf[0] && vc.isBiallelic())
			array.incr(VariantEnum.TRANSVERSIONS);
	}

	private boolean isTransition(Allele ref, Allele alt) throws IllegalArgumentException {
		byte refAllele = ref.getBases()[0];
		byte altAllele = alt.getBases()[0];
		return (refAllele == 'A' && altAllele == 'G') || (refAllele == 'G' && altAllele == 'A')
				|| (refAllele == 'C' && altAllele == 'T') || (refAllele == 'T' && altAllele == 'C');
	}

	private boolean[] basicStatics(Genotype gt, IntArray array, Allele Ref,GenotypeType type) {
		boolean[] conf = new boolean[4];
		Arrays.fill(conf, false);

		List<Allele> alts = gt.getAlleles();
		if(type == GenotypeType.NO_CALL || type == GenotypeType.HOM_REF)
			return conf;
		boolean isBiallelic = determinePolymorphicType(Ref, alts, conf);
		if(!isBiallelic)
			return new boolean[]{false,false,false,false};
		
		array.incr(VariantEnum.VARIANTS);
		if (conf[0])
			array.incr(VariantEnum.SNPS);
		if (conf[1])
			array.incr(VariantEnum.MNPS);
		if (conf[2]) {
			array.incr(VariantEnum.INSERTIONS);
		}
		if (conf[3]) {
			array.incr(VariantEnum.DELETIONS);
		}

		if (conf[2] || conf[3])
			array.incr(VariantEnum.INDELS);

		if (conf[0] && isBiallelic && isTransition(alts.get(0),alts.get(1))) {
			array.incr(VariantEnum.TRANSITIONS);
		} else if (conf[0] && isBiallelic)
			array.incr(VariantEnum.TRANSVERSIONS);

		return conf;
	}

	private boolean determinePolymorphicType(Allele Ref, List<Allele> alts, boolean[] conf) {
		boolean isBiallelic = alts.size() == 2 ? true : false;

		Allele ref = alts.get(1).isReference() ? alts.get(1) : alts.get(0);
		Allele allele = alts.get(1).isReference() ? alts.get(0) : alts.get(1);

		if (ref.isSymbolic() || allele.isSymbolic())
			return isBiallelic;

		if (ref.length() == allele.length()) {
			if (allele.length() == 1)
				conf[0] = true;
			else
				conf[1] = true;
		} else if (Ref.equals(ref)) {
			if (ref.length() > allele.length())
				conf[3] = true;
			else
				conf[2] = true;
		} else {
			if (ref.length() > Ref.length() || allele.length() > Ref.length())
				conf[2] = true;
			else
				conf[3] = true;
		}

		return isBiallelic;
	}
	
	private GenotypeType determineType(final List<Allele> alleles) {
        if ( alleles.isEmpty() )
            return GenotypeType.UNAVAILABLE;

        boolean sawNoCall = false, sawMultipleAlleles = false;
        Allele observedAllele = null;

        for ( final Allele allele : alleles ) {
            if ( allele.isNoCall() )
                sawNoCall = true;
            else if ( observedAllele == null )
                observedAllele = allele;
            else if ( !allele.equals(observedAllele) )
                sawMultipleAlleles = true;
        }

        if ( sawNoCall ) {
            if ( observedAllele == null )
                return GenotypeType.NO_CALL;
            return GenotypeType.MIXED;
        }

        if ( observedAllele == null )
            throw new IllegalStateException("BUG: there are no alleles present in this genotype but the alleles list is not null");

        return sawMultipleAlleles ? GenotypeType.HET : observedAllele.isReference() ? GenotypeType.HOM_REF : GenotypeType.HOM_VAR;
    }

	private void genotype(Genotype gt, String contig, IntArray statics, boolean[] conf,GenotypeType type) {
		if(type == GenotypeType.NO_CALL)
			return;
		boolean snp = conf[0], mnp = conf[1], insertion = conf[2], deletion = conf[3];
		contig = contig.toUpperCase();
		
		if (type ==  GenotypeType.HOM_VAR || type == GenotypeType.HOM_REF) {
			statics.incr(VariantEnum.HOM);
			if (snp)
				statics.incr(VariantEnum.SNP_HOM);
			if (mnp)
				statics.incr(VariantEnum.MNP_HOM);
			if (insertion) {
				statics.incr(VariantEnum.INSERTION_HOM);
			}
			if (deletion) {
				statics.incr(VariantEnum.DELETION_HOM);
			}

			if (insertion || deletion)
				statics.incr(VariantEnum.INDEL_HOM);

			if (contig.equals("CHRX") || contig.equals("X")) {
				statics.incr(VariantEnum.CHROM_X_HOM);
			}

			if (type == GenotypeType.HOM_REF)
				statics.incr(VariantEnum.HOM_RR);
			else
				statics.incr(VariantEnum.HOM_AA);
		}else{
			statics.incr(VariantEnum.HET);
			if (snp)
				statics.incr(VariantEnum.SNP_HET);
			if (mnp)
				statics.incr(VariantEnum.MNP_HET);
			if (insertion) {
				statics.incr(VariantEnum.INSERTION_HET);
			}
			if (deletion) {
				statics.incr(VariantEnum.DELETION_HET);
			}

			if (insertion || deletion)
				statics.incr(VariantEnum.INDEL_HET);

			if (contig.equals("CHRX") || contig.equals("X")) {
				statics.incr(VariantEnum.CHROM_X_HET);
			}

			if (gt.getAllele(0).isNonReference() && gt.getAllele(1).isNonReference())
				statics.incr(VariantEnum.HET_AA);
			else
				statics.incr(VariantEnum.HET_RA);
		}
	}

	private void genotyStatic(VariantContext vc, boolean[] conf) {
		boolean multiSample = false;
		if (vc.getSampleNames().size() > 1) {
			multiSample = true;
			sampleName = "statistics";
		}
		for (Genotype gt : vc.getGenotypesOrderedByName()) {
			String name = gt.getSampleName();
			if (!multiSample) {
				sampleName = name;
				sampleStatics.put(name, basicArray);
			}

			if (!sampleStatics.containsKey(name)) {
				IntArray temp = new IntArray(ALL_LENGTH, 0);
				sampleStatics.put(name, temp);
			}

			IntArray statics = sampleStatics.get(name);
			GenotypeType type = determineType(gt.getAlleles());
			if (multiSample) {
				boolean[] temp = basicStatics(gt, statics, vc.getReference(),type);
				genotype(gt, vc.getContig(), statics, temp,type);
			} else
				genotype(gt, vc.getContig(), statics, conf,type);
		}
	}

	private String getRateValue(int left, int right) {
		if (right == 0)
			return "NA";
		return String.format("%.3f", (float) left / right);
	}

	private String getFloat(IntArray array, VariantEnum left, VariantEnum right) {
		return getRateValue(array.get(left), array.get(right));
	}

	static class StringFormat {
		public String mode;
		public String key;
		public String value;

		public StringFormat(String mode, String key, String value) {
			this.mode = mode;
			this.key = key;
			this.value = value;
		}

		public StringFormat(String mode, String key, int value) {
			this.mode = mode;
			this.key = key;
			this.value = String.valueOf(value);
		}
	}

	private String format(ArrayList<StringFormat> list) {
		return format(null, list);
	}

	private String format(String name, ArrayList<StringFormat> list) {
		StringBuilder sb = new StringBuilder();

		if (name != null)
			sb.append(name + ":\n");
		for (StringFormat format : list) {
			sb.append(String.format(format.mode, format.key, format.value));
		}

		return sb.toString();
	}

	private String format(String f, int left, int right) {
		return String.format("%s(%d/%d)", f, left, right);
	}

	private void addBasicStatics(ArrayList<StringFormat> list, IntArray array) {
		list.add(new StringFormat(mode1, "variants:", array.get(VariantEnum.VARIANTS)));
		list.add(new StringFormat(mode1, "snps:", array.get(VariantEnum.SNPS)));
		list.add(new StringFormat(mode1, "mnps:", array.get(VariantEnum.MNPS)));
		list.add(new StringFormat(mode1, "insertions:", array.get(VariantEnum.INSERTIONS)));
		list.add(new StringFormat(mode1, "deletions:", array.get(VariantEnum.DELETIONS)));
		list.add(new StringFormat(mode1, "indels:", array.get(VariantEnum.INDELS)));
		list.add(new StringFormat(mode1, "snp Transitions/Transversions rate:",
				format(getFloat(array, VariantEnum.TRANSITIONS, VariantEnum.TRANSVERSIONS),
						array.get(VariantEnum.TRANSITIONS), array.get(VariantEnum.TRANSVERSIONS))));
	}

	private void addExtendStatics(ArrayList<StringFormat> list, IntArray array) {
		list.add(new StringFormat(mode1, "het/hom rate:", format(getFloat(array, VariantEnum.HET, VariantEnum.HOM),
				array.get(VariantEnum.HET), array.get(VariantEnum.HOM))));
		list.add(new StringFormat(mode1, "snp het/hom rate:",
				format(getFloat(array, VariantEnum.SNP_HET, VariantEnum.SNP_HOM), array.get(VariantEnum.SNP_HET),
						array.get(VariantEnum.SNP_HOM))));
		list.add(new StringFormat(mode1, "mnp het/hom rate:",
				format(getFloat(array, VariantEnum.MNP_HET, VariantEnum.MNP_HOM), array.get(VariantEnum.MNP_HET),
						array.get(VariantEnum.MNP_HOM))));
		list.add(new StringFormat(mode1, "insertion het/hom rate:",
				format(getFloat(array, VariantEnum.INSERTION_HET, VariantEnum.INSERTION_HOM),
						array.get(VariantEnum.INSERTION_HET), array.get(VariantEnum.INSERTION_HOM))));
		list.add(new StringFormat(mode1, "deletion het/hom rate:",
				format(getFloat(array, VariantEnum.DELETION_HET, VariantEnum.DELETION_HOM),
						array.get(VariantEnum.DELETION_HET), array.get(VariantEnum.DELETION_HOM))));
		list.add(new StringFormat(mode1, "indel het/hom rate:",
				format(getFloat(array, VariantEnum.INDEL_HET, VariantEnum.INDEL_HOM), array.get(VariantEnum.INDEL_HET),
						array.get(VariantEnum.INDEL_HOM))));
		list.add(new StringFormat(mode1, "insertion/deletion rate:",
				format(getFloat(array, VariantEnum.INSERTIONS, VariantEnum.DELETIONS),
						array.get(VariantEnum.INSERTIONS), array.get(VariantEnum.DELETIONS))));
		int snps = array.get(VariantEnum.SNPS) + array.get(VariantEnum.MNPS);
		list.add(new StringFormat(mode1, "indel/snp+mnp rate:",
				format(getRateValue(array.get(VariantEnum.INDELS), snps), array.get(VariantEnum.INDELS), snps)));
		list.add(new StringFormat(mode1, "chrx het/hom rate:",
				format(getFloat(array, VariantEnum.CHROM_X_HET, VariantEnum.CHROM_X_HOM),
						array.get(VariantEnum.CHROM_X_HET), array.get(VariantEnum.CHROM_X_HOM))));
		list.add(new StringFormat(mode1, "hom_AA/hom_RR rate:",
				format(getFloat(array, VariantEnum.HOM_AA, VariantEnum.HOM_RR), array.get(VariantEnum.HOM_AA),
						array.get(VariantEnum.HOM_RR))));
		list.add(new StringFormat(mode1, "het_AA/het_RA rate:",
				format(getFloat(array, VariantEnum.HET_AA, VariantEnum.HET_RA), array.get(VariantEnum.HET_AA),
						array.get(VariantEnum.HET_RA))));
	}

	public String toString() {
		ArrayList<StringFormat> list = new ArrayList<StringFormat>();
		StringBuilder sb = new StringBuilder();
		if(individuals && sampleStatics.size() > 1){
			addBasicStatics(list, basicArray);
			sb.append(format(list));
			list.clear();
		}
		if (individuals) {
			for (String name : sampleStatics.keySet()) {
				IntArray array = sampleStatics.get(name);
				addBasicStatics(list, array);
				if (array.length() > BASIC_LENGTH) {
					addExtendStatics(list, array);
				}
				sb.append(format(name, list));
				list.clear();
			}
		}
		return sb.toString();
	}

	public String getSampleName() {
		return this.sampleName;
	}
	
	public int getMixed(){
		return mixed;
	}/*

	public static void main(String[] args) throws IOException {
		VariantBasicStatistic statics = new VariantBasicStatistic(true);
		File f = new File("F:\\multi.vcf");
		AsciiLineReaderIterator iterator = null;
		if (f.getName().endsWith(".gz"))
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new GZIPInputStream(new FileInputStream(f))));
		else
			iterator = new AsciiLineReaderIterator(new AsciiLineReader(new FileInputStream(f)));
		VCFCodec codec = new VCFCodec();
		codec.readHeader(iterator);
		while (iterator.hasNext()) {
			VariantContext vc = codec.decode(iterator.next());
			statics.variantStatic(vc);
		}

		System.out.println(statics.toString());
	}*/
}
