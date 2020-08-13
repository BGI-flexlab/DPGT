package org.bgi.flexlab.gaea.tools.haplotypecaller.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.jointcalling.util.RodBinding;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.VariantContext;

public class RefMetaDataTracker {
	private GenomeLocation location = null;
	
	private HashMap<String,List<VariantContext>> dataSources = null;
	
	public static final String ALLELE_VALUE = "ALLELE_VALUE";
	public static final String DB_VALUE = "DB";
	
	public RefMetaDataTracker(GenomeLocation location){
		this.location = location;
	}
	
	public RefMetaDataTracker(){}
	
	
	public GenomeLocation getLocation(){
		return location;
	}
	
	public static List<VariantContext> getValues(List<VariantContext> binding,GenomeLocation location){
		if(location.getStart() != location.getStop())
			throw new UserException("location must length is 1!");
		List<VariantContext> results = new ArrayList<VariantContext>();
		int pos = location.getStart();
		for(VariantContext context : binding){
			if(context.getStart() == pos)
				results.add(context);
		}
		
		return results;
	}
	
	@SuppressWarnings("unchecked")
	public <T extends Feature> List<T> getValues(RodBinding<VariantContext> binding,GenomeLocation loc){
		return Collections.EMPTY_LIST;
	}
	
	public void add(String name,List<VariantContext> values){
		if(dataSources == null){
			dataSources = new HashMap<String,List<VariantContext>>();
		}
		
		if(dataSources.containsKey(name))
			throw new UserException("RefMetaDataTracker couldn't add same name="+name);
		
		dataSources.put(name, values);
	}
	
	public List<VariantContext> getValues(String name){
		if(dataSources == null || !dataSources.containsKey(name))
			return Collections.EMPTY_LIST;
		return dataSources.get(name);
	}
}
