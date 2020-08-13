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
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.structure.vcf;

import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import org.bgi.flexlab.gaea.util.ChromosomeUtils;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.zip.GZIPInputStream;

public abstract class AbstractVCFLoader {
	
	public static final Set<String> BLOCK_COMPRESSED_EXTENSIONS = Collections.unmodifiableSet(new HashSet<String>(Arrays.asList(".gz", ".gzip", ".bgz", ".bgzf")));
		
	protected String input;
			
	protected AsciiLineReaderIterator iterator;
	
	protected long pos;

	private VCFCodec codec;
	
	private FeatureCodecHeader header;
		
	public AbstractVCFLoader(String dbSNP) throws IOException {
		codec = new VCFCodec();
		input = dbSNP;
		readHeader();
		getSeekableStream();
		seek(0);
	}
	
	public abstract void getSeekableStream();
	
	public boolean hasNext() {
		return iterator.hasNext();
	}
	
	public PositionalVariantContext next() {
		long pos = iterator.getPosition() + this.pos;
		String line = iterator.next();
		while(line.startsWith(VCFHeader.HEADER_INDICATOR)) {
			pos = iterator.getPosition() + this.pos;
			line = iterator.next();
		}
		VariantContext vc = codec.decode(line);
		PositionalVariantContext pvc = new PositionalVariantContext(vc, pos);
		return pvc;
	}
	
	@Deprecated
	public Iterator<PositionalVariantContext> iterator() throws IOException {
		return collect();
	}
	
	@Deprecated
	public Iterator<PositionalVariantContext> iterator(long pos) throws IOException {
		seek(pos);;
		return collect();
	}
	
	private Iterator<PositionalVariantContext> collect() throws IOException {
		ArrayList<PositionalVariantContext> variantContexts = new ArrayList<>();
		while(iterator.hasNext()) {
			long pos = iterator.getPosition() + this.pos;
			String line = iterator.next();
			if(line.startsWith(VCFHeader.HEADER_INDICATOR))
				continue;
			VariantContext vc = codec.decode(line);
			variantContexts.add(new PositionalVariantContext(vc, pos));
		}
		return variantContexts.iterator();
	}
	
	public static class PositionalVariantContext {
		
		private VariantContext vc;
		
		private long position;
		
		public PositionalVariantContext(VariantContext vc, long position) {
			this.vc = vc;
			this.position = position;
		}
		
		public VariantContext getVariantContext() {
			return vc;
		}
		
		public long getPosition() {
			return position;
		}
		
	}
	
	
	public ArrayList<VariantContext> query(String chrName, long start, int end) throws IOException {
		ArrayList<VariantContext> contexts = new ArrayList<VariantContext>();

		chrName = ChromosomeUtils.formatChrName(chrName);
		seek(start);

		while (hasNext()) {
			VariantContext context = next().getVariantContext();
			if (chrName.equals(ChromosomeUtils.formatChrName(context.getContig()))) {
				if (context.getStart() > end) {
					break;
				}
				contexts.add(context);
			} else {
				break;
			}
		}

		return contexts;
	}

//	code derived from htsjdk
	private void readHeader() throws IOException {
		InputStream is = null;
        PositionalBufferedStream pbs = null;
        try {
            is = getInputStream(input);
            if (hasBlockCompressedExtension(input)) {
                // TODO -- warning I don't think this can work, the buffered input stream screws up position
                is = new GZIPInputStream(new BufferedInputStream(is));
            }
            pbs = new PositionalBufferedStream(is);
            header = codec.readHeader(codec.makeSourceFromStream(pbs));
        } catch (Exception e) {
            throw new TribbleException.MalformedFeatureFile("Unable to parse header with error: " + e.getMessage(), input, e);
        } finally {
            if (pbs != null) pbs.close();
            else if (is != null) is.close();
        }	
	}
	
	public abstract InputStream getInputStream (String input) throws IOException;
//	public static String format(String inputFile) {
//		return inputFile + VcfIndex.INDEX_SUFFIX;
//	}
	
	public VCFHeader getHeader() {
		return (VCFHeader)(header.getHeaderValue());
	}
	
	public abstract void seek(long pos);
	
	public static boolean hasBlockCompressedExtension(String fileName){
		 for (final String extension : BLOCK_COMPRESSED_EXTENSIONS) {
	            if (fileName.toLowerCase().endsWith(extension))
	                return true;
	     }
	     return false;
	}
	
	public void close() {
		if(iterator != null)
			try {
				iterator.close();
			} catch (IOException e) {
				throw new RuntimeException(e.toString());
			}
	}
	
	public static boolean isVCFStream(final InputStream stream, final String MAGIC_HEADER_LINE) {
        try {
            byte[] buff = new byte[MAGIC_HEADER_LINE.length()];
            int nread = stream.read(buff, 0, MAGIC_HEADER_LINE.length());
            boolean eq = Arrays.equals(buff, MAGIC_HEADER_LINE.getBytes());
            return eq;
        } catch ( IOException e ) {
            return false;
        } catch ( RuntimeException e ) {
            return false;
        } finally {
            try { stream.close(); } catch ( IOException e ) {}
        }
    }
	
}
