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
package org.bgi.flexlab.gaea.tools.vcfqualitycontrol;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class HardFilter {
	private List<JexlVCMatchExp> snpFilter = null;
	private List<JexlVCMatchExp> indelFilter = null;
	protected Boolean FAIL_MISSING_VALUES = false;
	public HardFilter(String filterName, String snpFilter, String indelFilter) {
		if(snpFilter != null)
			this.snpFilter = VariantContextUtils.initializeMatchExps(new String[]{filterName+"snp"}, new String[]{snpFilter});
		if(indelFilter != null)
			this.indelFilter = VariantContextUtils.initializeMatchExps(new String[]{filterName+"indel"}, new String[]{indelFilter});
		VariantContextUtils.engine.get().setSilent(true);
	}
	
	public VCFHeader addFilterHeader(VCFHeader header) {
		for(VCFHeaderLine fLine : filterHeaders())
			header.addMetaDataLine(fLine);
		return header;
	}
	
	public Set<VCFHeaderLine> filterHeaders() {
		Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
		for ( VariantContextUtils.JexlVCMatchExp exp : snpFilter )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
		for ( VariantContextUtils.JexlVCMatchExp exp : indelFilter )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
		return hInfo;
	}
	
	public VariantContext filter(VariantContext vc) {
		final VariantContextBuilder builder = new VariantContextBuilder(vc);
		Set<String> filters = new LinkedHashSet<String>(vc.getFilters());
		if(vc.isSNP() && snpFilter != null) {
			for ( VariantContextUtils.JexlVCMatchExp exp : snpFilter ) {
				try {
	                if ( VariantContextUtils.match(vc, exp) ) {
	                	filters.add(exp.name);
	                }
	            } catch (Exception e) {
	                // do nothing unless specifically asked to; it just means that the expression isn't defined for this context
	                e.printStackTrace();
	            	if ( FAIL_MISSING_VALUES )
	                    filters.add(exp.name);                         
	            }
			}
		}
		
		if(vc.isIndel() && indelFilter != null) {
			for ( VariantContextUtils.JexlVCMatchExp exp : indelFilter ) {
				try {
	                if ( VariantContextUtils.match(vc, exp) )
	                    filters.add(exp.name);
	            } catch (Exception e) {
	                // do nothing unless specifically asked to; it just means that the expression isn't defined for this context
	                if ( FAIL_MISSING_VALUES )
	                    filters.add(exp.name);                         
	            }
			}
		}
		if ( filters.isEmpty() )
            builder.passFilters();
        else 
            builder.filters(filters);
        
		return builder.make();
	}
}
