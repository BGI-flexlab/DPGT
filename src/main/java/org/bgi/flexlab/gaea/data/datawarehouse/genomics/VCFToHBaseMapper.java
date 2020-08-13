package org.bgi.flexlab.gaea.data.datawarehouse.genomics;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.client.Put;
import org.apache.hadoop.hbase.io.ImmutableBytesWritable;
import org.apache.hadoop.hbase.util.Bytes;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.header.GaeaVCFHeader;
import org.bgi.flexlab.gaea.data.structure.header.MultipleVCFHeader;
import org.bgi.flexlab.gaea.util.Utils;
import org.seqdoop.hadoop_bam.VariantContextWritable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
public class VCFToHBaseMapper extends Mapper<LongWritable, VariantContextWritable, ImmutableBytesWritable, Put> {
	private Set<String> sampleNames = new TreeSet<String>();

	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		Configuration conf = context.getConfiguration();
		MultipleVCFHeader headers = (MultipleVCFHeader) GaeaVCFHeader.loadVcfHeader(false, conf);

		for (String sample : headers.getSampleNames(0)) {
			sampleNames.add(sample);
		}
	}
	
	public HashMap<Allele, Integer> buildAlleleStrings(final VariantContext vc) {
		final HashMap<Allele, Integer> alleleMap = new HashMap<Allele, Integer>(vc.getAlleles().size()+1);
		alleleMap.put(Allele.NO_CALL, 3);

		final List<Allele> alleles = vc.getAlleles();
		for ( int i = 0; i < alleles.size(); i++ ) {
			alleleMap.put(alleles.get(i), i);
		}

		return alleleMap;
	}
	
	private void createRecord(String rowKey,String ref,String alt,String quality,byte gt,String ad,int dp,Context context ) throws IOException, InterruptedException{
		byte[] rowKeys = Bytes.toBytes(rowKey);
		
		byte[] gts = new byte[1];
		gts[0] = gt;
		
		context.getCounter("Counters", "record").increment(1);
		
		ImmutableBytesWritable rKey = new ImmutableBytesWritable(rowKeys);
		Put put = new Put(rowKeys);
		put.add(Bytes.toBytes("info"),Bytes.toBytes("REF"),Bytes.toBytes(ref));
		put.add(Bytes.toBytes("info"),Bytes.toBytes("ALT"),Bytes.toBytes(alt));
		put.add(Bytes.toBytes("info"),Bytes.toBytes("QUAL"),Bytes.toBytes(quality));
		put.add(Bytes.toBytes("info"),Bytes.toBytes("GT"),gts);
		if(ad != null && ad.trim().length() > 0)
			put.add(Bytes.toBytes("info"),Bytes.toBytes("AD"),Bytes.toBytes(ad));
		if(dp >= 0)
			put.add(Bytes.toBytes("info"),Bytes.toBytes("DP"),Bytes.toBytes(dp));
		context.write(rKey, put);
	}
	
	private String formatQualValue(final double qual) {
		String s = String.format("%.2f", qual);
		if ( s.endsWith(".00") )
			s = s.substring(0, s.length() - ".00".length());
		return s;
	}
	
	private void parserVariantContext(VariantContext vc,Context context,HashMap<Allele,Integer> alleleMap) throws IOException, InterruptedException{
		HashMap<String,Integer> altMap = new HashMap<String,Integer>();
		ArrayList<String> list = new ArrayList<String>();
		final int ploidy = vc.getMaxPloidy(2);
		int index;
		int currentIndex = 1;
		boolean nonHomRef = false;
		
		String rowKey = vc.getContig()+"-"+vc.getStart()+"-";
		String quality = formatQualValue(vc.getPhredScaledQual());
		String ref = vc.getReference().getDisplayString();

		for (final String sample : sampleNames) {
			currentIndex = 1;
			altMap.clear();
			list.clear();
			nonHomRef = false;
			Genotype g = vc.getGenotype(sample);
			if (g == null) 
				g = GenotypeBuilder.createMissing(sample, ploidy);
			
			// g.getPloidy() must be equals to two
			for (int i = 0; i < g.getPloidy(); i++) {
				if(!alleleMap.containsKey(g.getAllele(i)))
					throw new UserException("Variant conext alleles is not contains genotype allele");
				index = alleleMap.get(g.getAllele(i));
				if(index == 0)
					altMap.put(g.getAllele(i).getDisplayString(), index);
				else{
					nonHomRef = true;
					if(!altMap.containsKey(g.getAllele(i).getDisplayString())){
						altMap.put(g.getAllele(i).getDisplayString(), currentIndex);
						list.add(g.getAllele(i).getDisplayString());
						currentIndex++;
					}
				}
			}
			
			if(nonHomRef){
				String alt = Utils.join(",", list);
				if(alt.trim().length() == 0)
					continue;
				int gt = 0;
				for(int i = 0 ; i < g.getPloidy();i++){
					gt <<= 4;
					gt |= (altMap.get(g.getAllele(i).getDisplayString()) & 0xf);
				}
				
				createRecord(rowKey+sample,ref,alt,quality,(byte)(gt&0xff),Utils.join(",", g.getAD()),g.getDP(),context);
			}
		}
		
		altMap.clear();
		list.clear();
	}

	@Override
	public void map(LongWritable key, VariantContextWritable value, Context context)
			throws IOException, InterruptedException {
		VariantContext ctx = value.get();
		
		if(ctx.isVariant()){
			parserVariantContext(ctx,context,buildAlleleStrings(ctx));
		}
	}
}
