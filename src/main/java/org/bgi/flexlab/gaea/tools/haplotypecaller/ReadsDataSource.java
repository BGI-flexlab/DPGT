package org.bgi.flexlab.gaea.tools.haplotypecaller;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pileup.ReadsDownsampler;
import org.bgi.flexlab.gaea.tools.haplotypecaller.readfilter.ReadFilter;

import htsjdk.samtools.SAMFileHeader;

public class ReadsDataSource implements Iterator<GaeaSamRecord>, Iterable<GaeaSamRecord> {
	private GaeaSamRecord currentRecord = null;

	// reads iterator
	private Iterator<SamRecordWritable> reads;

	private List<GaeaSamRecord> overlaps = new ArrayList<GaeaSamRecord>();

	private SAMFileHeader header = null;

	private ReadsDownsampler downsampler;
	private Iterator<GaeaSamRecord> cachedDownsampledReads = null;
	private GaeaSamRecord nextRead = null;
	private GenomeLocation queryInterval = null;
	private ReadFilter readFilter = null;
	
	private int readsNumber = 0;

	public ReadsDataSource(Iterable<SamRecordWritable> iterable, SAMFileHeader header) {
		this.header = header;
		dataReset(iterable);
	}

	public void dataReset(Iterable<SamRecordWritable> iterable) {
		reads = iterable.iterator();
		if (overlaps != null && !overlaps.isEmpty())
			overlaps.clear();
		init();
	}

	private void init() {
		if (reads.hasNext())
			currentRecord = new GaeaSamRecord(header, reads.next().get());
		else
			currentRecord = null;
		
		if(currentRecord != null)
			readsNumber = 1;
		readsNumber = 0;
	}

	private int overlapReads(GaeaSamRecord read) {
		if (read.getEnd() < queryInterval.getStart()) {
			return 0;
		} else if ((read.getAlignmentStart() >= queryInterval.getStart()
				&& read.getAlignmentStart() <= queryInterval.getEnd())
				|| (read.getAlignmentEnd() >= queryInterval.getStart()
						&& read.getAlignmentEnd() <= queryInterval.getEnd())) {
			return 1;
		} else
			return 2;
	}

	private void processReads() {
		final boolean traversalIsBounded = queryInterval != null;
		if (traversalIsBounded) {
			while (currentRecord != null) {
				if ((readFilter != null && readFilter.test(currentRecord)) || readFilter == null) {
					int result = overlapReads(currentRecord);
					if (result == 1) {
						if (downsampler != null) {
							downsampler.submit(currentRecord);
							
						}
						else
							overlaps.add(currentRecord);
					} else if (result == 2)
						break;
				}
				if (reads.hasNext()) {
					currentRecord = new GaeaSamRecord(header, reads.next().get());
					readsNumber++;
				}
				else
					currentRecord = null;
			}
		}
	}
	
	private void processPrevReads() {
		if (overlaps != null && !overlaps.isEmpty()) {
			ArrayList<GaeaSamRecord> removes = new ArrayList<GaeaSamRecord>();
			for (GaeaSamRecord read : overlaps) {
				if (overlapReads(read) == 1) {
					if (downsampler != null)
						downsampler.submit(read);
				} else
					removes.add(read);
			}

			if (downsampler != null)
				overlaps.clear();
			else {
				overlaps.removeAll(removes);
				removes.clear();
			}
		}
	}

	private boolean fillDownsampledReadsCache() {
		processPrevReads();
		processReads();

		if (downsampler != null) {
			//if (!reads.hasNext())
				downsampler.signalEndOfInput();
			overlaps = downsampler.consumeFinalizedItems();
			downsampler.clearItems();
		}
		cachedDownsampledReads = overlaps.iterator();

		return cachedDownsampledReads.hasNext();
	}

	public void clear() {
		if (overlaps != null)
			overlaps.clear();
		if(downsampler != null)
			downsampler.clearItems();
		nextRead = null;
	}

	@Override
	public Iterator<GaeaSamRecord> iterator() {
		return this;
	}

	@Override
	public boolean hasNext() {
		return nextRead != null;
	}

	@Override
	public GaeaSamRecord next() {
		if (nextRead == null) {
			throw new NoSuchElementException("next() called when there are no more items");
		}

		final GaeaSamRecord toReturn = nextRead;
		nextReads();

		return toReturn;
	}

	private void nextReads() {
		if (readyToReleaseReads())
			nextRead = cachedDownsampledReads.next();
		else
			nextRead = null;
	}

	private void advanceToNextRead() {
		if (fillDownsampledReadsCache()) {
			nextRead = cachedDownsampledReads.next();
		} else {
			nextRead = null;
		}
	}

	public void set(GenomeLocation interval, ReadFilter readFilter, ReadsDownsampler downsampler) {
		this.queryInterval = interval;
		this.readFilter = readFilter;
		this.downsampler = downsampler;
		cachedDownsampledReads = null;
		advanceToNextRead();
	}

	private boolean readyToReleaseReads() {
		return cachedDownsampledReads != null && cachedDownsampledReads.hasNext();
	}
	
	public void printReads() {
		if(queryInterval.getContig().equals("chr1") && queryInterval.getStart() == 12600 && queryInterval.getEnd() == 13099) {
			for(GaeaSamRecord read : overlaps) {
				System.err.println(read.toString());
			}
		}
	}
	
	public int getReadsNumber() {
		return this.readsNumber;
	}
}
