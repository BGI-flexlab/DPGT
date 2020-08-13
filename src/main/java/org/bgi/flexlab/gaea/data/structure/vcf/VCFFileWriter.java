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

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.hadoop.conf.Configuration;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;

public abstract class VCFFileWriter implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -2686163989598391782L;
	
	private VariantContextWriterBuilder builder;
	
	private VariantContextWriter writer;
		
	private boolean hasWrittenHeader = false;
	
	protected OutputStream os;
	
	public VCFFileWriter(String path,boolean doNotWriteGenotypes,
			final boolean allowMissingFieldsInHeader, Configuration conf) throws IOException {
		initOutputStream(path, conf);
		initBuilder(os, doNotWriteGenotypes, allowMissingFieldsInHeader);
		this.writer = builder.build();
	}
	
	public abstract void initOutputStream(String path, Configuration conf); 
	
	private void initBuilder(OutputStream os, boolean doNotWriteGenotypes, 
			boolean allowMissingFieldsInHeader) {
		builder = new VariantContextWriterBuilder();
		builder.clearOptions();
		if(doNotWriteGenotypes) {
			builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);
		}
		if(allowMissingFieldsInHeader) {
			builder.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
		}
		builder.setOutputStream(os);
	}
	
	public void writeHeader(VCFHeader header) {
		writer.writeHeader(header);
		hasWrittenHeader = true;
	}
	
	public void close() {
		writer.close();
	}
	
	public void add(VariantContext vc) {
		if (!hasWrittenHeader) {
			throw new IllegalStateException(
					"The VCF Header must be written before records can be added: ");
		}
		writer.add(vc);
	}
}
