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
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.vcf;

import htsjdk.tribble.index.Index;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.bgi.flexlab.gaea.data.mapreduce.input.vcf.VCFHdfsLoader;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.vcf.index.IndexCreator;
import org.bgi.flexlab.gaea.util.Window;

import java.io.IOException;
import java.util.*;



public class VariantDataTracker {
	
	public enum Style {
		ALL,
		POS;
	}
	
	private String name;
	private String type;
	private Window window;
	private VCFHdfsLoader loader = null;;
	private List<VariantContext> data = null;
	private HashSet<Integer> site = null;
	private boolean bound;
	private Style style;
	
	public VariantDataTracker(){
		this.bound=false;
	};
	
	public VariantDataTracker(String source,String type ,Style s) {
		this.name=source;
		this.type=type;
		this.bound=true;
		this.style=s;
		initializeLoader();
	}
	
	public VariantDataTracker(String source,String type, Window window, Style s) {
		this(source, type, s);
		this.window=window;
		initializeData();
	}
	
	public VariantDataTracker(VCFHdfsLoader loader,String type,Style s) {
		this.loader=loader;
		this.type=type;
		this.style=s;
	}
	
	public void setSource(String source,String type) throws IOException {
		this.name=source;
		this.type=type;
		initializeLoader();
		initializeData();
	}
	
	/*public void setWindow(Window win) throws IOException
	{
		this.window=win;
		initializeData();
	}*/
	
	private void initializeLoader() {
		if(loader==null) {
			//creator.initialize(name, creator.defaultBinSize());
			IndexCreator creator = new IndexCreator(name);
			try {
				Index idx = creator.finalizeIndex();
				loader = new VCFHdfsLoader(name);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}	
	}

	public void initializeData() {
		if(window == null)
			return;
		if(data == null)
			data = new ArrayList<VariantContext>();
		
		loadVariantContext();
		loadVariantContextPos();
	}

	public boolean getValue(GenomeLocation onlyAtThisLoc) {
		return addValue(onlyAtThisLoc, true, false);
	}

	private boolean addValue(GenomeLocation curLocation,  boolean requireStartHere,
			final boolean takeFirstOnly)  {
		HashSet<Integer> val=null;
		if(site!=null) {
			val = site;
		} else {
			try {
				List<VariantContext> temp = loader.query(curLocation.getContig(), curLocation.getStart(), curLocation.getStop());
				if(temp!=null) {
					val=new HashSet<Integer>();
					for(VariantContext c:temp) {
						val.add(c.getStart());
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		if(val==null || val.size() == 0)
			return false;
		if(!requireStartHere || val.contains(curLocation.getStart()))
			return true;
		return false;
	}
	
	public List<VariantContext> getValues(GenomeLocation onlyAtThisLoc) {
		return addValues(new ArrayList<VariantContext>(1), onlyAtThisLoc, true, false);
	}

	private List<VariantContext> addValues( List<VariantContext> values, GenomeLocation curLocation,  boolean requireStartHere,
			final boolean takeFirstOnly)  {
		List<VariantContext> val=null;
		val = loadVariantContext(val, curLocation);
		if(val==null)
			return Collections.<VariantContext> emptyList();
		
		for (VariantContext rec : val) {
			if (!requireStartHere || rec.getStart() == curLocation.getStart()) { 
				if (takeFirstOnly) {
					if (values == null)
						values = Arrays.asList(rec);
					else
						values.add(rec);
					break;
				} else {
					if (values == null)
						values = new ArrayList<VariantContext>();
					values.add(rec);
				}
			}
		}

		return values == null ? Collections.<VariantContext> emptyList() : values;
	}

	public VCFHeader getHeader() {
		if(loader!=null)
			return loader.getHeader();
		return null;
	}
	
	public String getName() {
		return name;
	}

	public String getType() {
		return type;
	}

	public boolean isBound() {
		return bound;
	}
	
	
	private void loadVariantContext() {
		try {
			ArrayList<VariantContext> temp = loader.query(window.getContigName(), window.getStart(), window.getStop());
			if(temp!=null)
				data = temp;
			else
				System.out.println("temp is null");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private List<VariantContext> loadVariantContext(List<VariantContext> val, GenomeLocation curLocation) {
		List<VariantContext> result = val;
		if(data!=null) {
			result = data;
		} else {
			try {
				result = loader.query(curLocation.getContig(), curLocation.getStart(), curLocation.getStop());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return result;
	}
	
	private void loadVariantContextPos() {
		if(style == Style.POS) {
			site = new HashSet<Integer>();
			if(data==null)
				return;
			
			for(VariantContext c:data) {
				site.add(c.getStart());
			}
			data=null;
		}
	}
}
