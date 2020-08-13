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
package org.bgi.flexlab.gaea.data.structure.vcf.report;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;

import java.util.HashSet;
import java.util.Set;

public class ReportDatum {
	private int snpCount;
	private int indelCount;
	private int ti;
	private int tv;
	private String snpChr;
	private static Set<String> snpChrs = new HashSet<>();
	private String indelChr;
	private static Set<String> indelChrs = new HashSet<>();
	private ReportDatum(Builder builder) {
		this.snpCount = builder.snpCount;
		this.indelCount = builder.indelCount;
		this.ti = builder.ti;
		this.tv = builder.tv;
		this.snpChr = builder.snpChr;
		this.indelChr = builder.indelChr;
		snpChrs.add(this.snpChr);
		indelChrs.add(this.indelChr);
	}
	
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append(snpCount);builder.append("\t");
		builder.append(ti);builder.append("\t");
		builder.append(tv);builder.append("\t");
		builder.append(indelCount);builder.append("\t");
		builder.append(snpChr);builder.append("\t");
		builder.append(indelChr);builder.append("\t");
		return builder.toString();
	}
	
	public String formatReport() {
		StringBuilder tbBuilder = new StringBuilder();
		tbBuilder.append("SNP counts: " + snpCount);tbBuilder.append("\n");
		tbBuilder.append("Indel counts: " + indelCount);tbBuilder.append("\n");
		tbBuilder.append("Ti/Tv: " + (ti * 1.0 / tv));tbBuilder.append("\n");
		tbBuilder.append("SNP chromosome counts: " + snpChrs.size());tbBuilder.append("\n");
		tbBuilder.append("Indel chromosome counts: " + indelChrs.size());tbBuilder.append("\n");
		return tbBuilder.toString();
	}
	
	public ReportDatum combine(ReportDatum vi) {
		this.snpCount += vi.getSnpCount();
		this.indelCount += vi.getIndelCount();
		this.ti += vi.getTi();
		this.tv += vi.getTv();
		return this;
	}
	
	public int getSnpCount() {
		return snpCount;
	}
	
	public int getIndelCount() {
		return indelCount;
	}
	
	public int getTi() {
		return ti;
	}
	
	public int getTv() {
		return tv;
	}
	
	public String getSnpChr() {
		return snpChr;
	}
	
	public String getIndelChr() {
		return indelChr;
	}
	
	public static class Builder {
		
		private int snpCount;
		private int indelCount;
		private int ti;
		private int tv;
		private String snpChr;
		private String indelChr;
		private VariantContext vc;
		
		public Builder(VariantContext vc) {
			this.vc = vc;
		}
		
		public Builder() {
			
		}
		
		public Builder isSnp() {
			if(vc.isSNP() && vc.isBiallelic()) {
				snpCount++;
				snpChr = vc.getChr();
			}
			return this;
		}
		
		public Builder isTransition() {
			if((snpCount != 0) && VariantContextUtils.isTransition(vc))
				ti++;
			else
				tv++;
			return this;
		}
		
		public Builder isIndel() {
			if(vc.isIndel())  {
				indelCount++;
				indelChr = vc.getChr();
			}
			return this;
		}
		
		public ReportDatum buildFrom(String context) {
			String[] info = context.split("\t");
			snpCount = Integer.parseInt(info[0]);
			ti = Integer.parseInt(info[1]);
			tv = Integer.parseInt(info[2]);
			indelCount = Integer.parseInt(info[3]);
			snpChr = info[4];
			indelChr = info[5];
			return new ReportDatum(this);
		}
		
		public ReportDatum build() {
			return new ReportDatum(this);
		}
	}
}
