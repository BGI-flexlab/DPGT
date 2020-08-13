package org.bgi.flexlab.gaea.data.structure.variant;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextTracker {
	private VariantContext curr_variant = null;
	private HashMap<String, VariantContextBinding> bindings = null;
	private Iterator iterator = null;

	public VariantContextTracker(Iterator iterator) {
		this.iterator = iterator;
		bindings = new HashMap<String, VariantContextBinding>();
	}

	public void clear() {
		if (bindings != null)
			bindings.clear();
		curr_variant = null;
	}
	
	public void lazyLoad(GenomeLocation location){
		if (curr_variant == null) {
			if (iterator.hasNext())
				curr_variant = (VariantContext) iterator.next();
			else
				return;
		}

		while (curr_variant != null && curr_variant.getEnd() < location.getStart()) {
			if (iterator.hasNext())
				curr_variant = (VariantContext) iterator.next();
			else
				curr_variant = null;
		}
		
		if(curr_variant == null)
			return;
		
		while(curr_variant != null && curr_variant.getStart() <= location.getStop()){
			String name = curr_variant.getID();
			
			VariantContextBinding binding = null;
			if(bindings.containsKey(name)){
				binding = bindings.get(name);
			}else{
				binding = new VariantContextBinding();
			}
			
			binding.add(curr_variant);
			
			if (iterator.hasNext())
				curr_variant = (VariantContext) iterator.next();
			else
				curr_variant = null;
		}
	}
	
	public List<VariantContext> getPrioritizedValue(GenomeLocation interval,boolean firstOnly){
		List<VariantContext> intervals = new ArrayList<VariantContext>();
		
		for(String sample : bindings.keySet()){
			VariantContextBinding binding = bindings.get(sample);
			if(firstOnly){
				VariantContext value = binding.getFirstValue(interval);
				if(value != null)
					intervals.add(null);
			}else{
				List<VariantContext> values = binding.getPrioritizedValue(interval);
				if(values != null && values.size() > 0)
					intervals.addAll(values);
			}
		}
		
		if(intervals.size() == 0)
			return null;
		
		return intervals;
	}

	public int depth(GenomeLocation location) {
		lazyLoad(location);
		
		int depth = 0;
		for(String sample : bindings.keySet()){
			VariantContextBinding bind = bindings.get(sample);
			bind.exclusion(location);
			depth += bind.size();
		}
		
		return depth;
	}
}
