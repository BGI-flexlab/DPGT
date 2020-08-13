package org.bgi.flexlab.gaea.data.structure.variant;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextBinding {
	private List<VariantContext> binding = null;
	
	public VariantContextBinding(){
		binding = new ArrayList<VariantContext>();
	}
	
	public void add(VariantContext ctx){
		binding.add(ctx);
	}
	
	public List<VariantContext> getPrioritizedValue(GenomeLocation interval){
		List<VariantContext> intervals = new ArrayList<VariantContext>();
		
		int size = binding.size();
		for(int i = 0 ; i < size ; i++){
			VariantContext ctx = binding.get(i);
			
			if(ctx.getStart() > interval.getStop())
				break;
			if(ctx.getEnd() < interval.getStart())
				continue;
			
			intervals.add(ctx);
		}
		
		if(intervals.size() == 0)
			return null;
		
		return intervals;
	}
	
	public VariantContext getFirstValue(GenomeLocation interval){
		VariantContext result = null;
		
		int size = binding.size();
		for(int i = 0 ; i < size ; i++){
			VariantContext ctx = binding.get(i);
			
			if(ctx.getStart() > interval.getStop())
				break;
			if(ctx.getEnd() < interval.getStart())
				continue;
			
			result = ctx;
			break;
		}
		
		return result;
	}
	
	public void exclusion(GenomeLocation interval){
		if(binding == null)
			return;
		
		Iterator<VariantContext> iterator = binding.iterator();
		
		while(iterator.hasNext()){
			VariantContext ctx = iterator.next();
			
			if(ctx.getEnd() < interval.getStart())
				iterator.remove();
			else
				break;
		}
	}
	
	public int size(){
		return binding.size();
	}
	
	public void clear(){
		if(binding != null)
			binding.clear();
	}
}
