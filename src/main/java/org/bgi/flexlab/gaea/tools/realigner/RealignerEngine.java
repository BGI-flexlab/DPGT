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
package org.bgi.flexlab.gaea.tools.realigner;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.dbsnp.DbsnpShare;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.reference.ChromosomeInformationShare;
import org.bgi.flexlab.gaea.data.structure.reference.ReferenceShare;
import org.bgi.flexlab.gaea.data.structure.reference.index.VcfIndex;
import org.bgi.flexlab.gaea.data.structure.vcf.VCFLocalLoader;
import org.bgi.flexlab.gaea.data.variant.filter.VariantRegionFilter;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RealignerOptions;
import org.bgi.flexlab.gaea.util.Window;

import java.util.ArrayList;

public class RealignerEngine {
	private RealignerOptions option = null;
	private ReferenceShare genomeShare = null;
	private VCFLocalLoader loader = null;
	private ChromosomeInformationShare chrInfo = null;
	private ArrayList<VariantContext> knowIndels = null;
	private ArrayList<GaeaSamRecord> records = null;
	private ArrayList<GaeaSamRecord> filterRecords = null;
	private Window win = null;
	private SAMFileHeader mHeader = null;
	private VariantRegionFilter indelFilter = null;
	private IndelRealigner indelRealigner = null;
	private RealignerWriter writer = null;
	private DbsnpShare dbsnpShare = null;

	private int start;
	private int end;

	public RealignerEngine(RealignerOptions option, ReferenceShare genomeShare, DbsnpShare dbsnpShare,
			VCFLocalLoader loader, SAMFileHeader mHeader, RealignerWriter writer) {
		this.option = option;
		this.genomeShare = genomeShare;
		this.loader = loader;
		this.mHeader = mHeader;
		this.writer = writer;
		this.dbsnpShare = dbsnpShare;
	}

	public void set(Window win, ArrayList<GaeaSamRecord> records, ArrayList<GaeaSamRecord> filterRecords) {
		this.win = win;
		if (win == null)
			throw new RuntimeException("window is null");
		this.records = records;
		this.filterRecords = filterRecords;
		indelFilter = new VariantRegionFilter();
		setChromosome(genomeShare);
		setKnowIndels(loader);
		indelRealigner = new IndelRealigner(mHeader, knowIndels, win, chrInfo, option);
	}

	private void setChromosome(ReferenceShare genomeShare) {
		chrInfo = genomeShare.getChromosomeInfo(win.getContigName());
	}

	private void setKnowIndels(VCFLocalLoader loader) {
		if (loader == null) {
			throw new RuntimeException("loader is null!!");
		}
		String referenceName = win.getContigName();
		int WINDOWS_EXTEND = option.getExtendSize();

		start = (win.getStart() - WINDOWS_EXTEND) > 0 ? (win.getStart() - WINDOWS_EXTEND) : 0;
		end = (win.getStop() + WINDOWS_EXTEND) < mHeader.getSequence(referenceName).getSequenceLength()
				? (win.getStop() + WINDOWS_EXTEND) : mHeader.getSequence(referenceName).getSequenceLength();

		long startPosition = dbsnpShare.getStartPosition(referenceName, start / VcfIndex.WINDOW_SIZE, end / VcfIndex.WINDOW_SIZE,
				VcfIndex.WINDOW_SIZE);

		if (startPosition >= 0)
			knowIndels = indelFilter.loadFilter(loader, referenceName, startPosition, end);
	}

	public void reduce() {
		IdentifyRegionsCreator creator = new IdentifyRegionsCreator(option, filterRecords, mHeader, chrInfo,
				knowIndels);
		creator.regionCreator(win.getChrIndex(), 0, Integer.MAX_VALUE);

		ArrayList<GenomeLocation> intervals = creator.getIntervals();
		filterRecords.clear();

		indelRealigner.setIntervals(intervals);
		indelRealigner.traversals(records, writer);
	}
}
