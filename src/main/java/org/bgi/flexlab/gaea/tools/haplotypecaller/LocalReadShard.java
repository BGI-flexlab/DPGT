package org.bgi.flexlab.gaea.tools.haplotypecaller;

import java.util.Iterator;
import java.util.List;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.tools.haplotypecaller.pileup.ReadsDownsampler;
import org.bgi.flexlab.gaea.tools.haplotypecaller.readfilter.ReadFilter;
import org.bgi.flexlab.gaea.util.Utils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;

public final class LocalReadShard implements Shard<GaeaSamRecord> {

	private final GenomeLocation interval;
	private final GenomeLocation paddedInterval;
	private final ReadsDataSource readsSource;
	private ReadFilter readFilter;
	private ReadsDownsampler downsampler;

	/**
	 * Create a new Shard spanning the specified interval, with the specified
	 * amount of padding.
	 *
	 * @param interval
	 *            the genomic span covered by this shard
	 * @param paddedInterval
	 *            the span covered by this shard, plus any additional padding on
	 *            each side (must contain the un-padded interval)
	 * @param readsSource
	 *            source of reads from which to populate this shard
	 */
	public LocalReadShard(final GenomeLocation interval, final GenomeLocation paddedInterval,
			final ReadsDataSource readsSource) {
		Utils.nonNull(interval);
		Utils.nonNull(paddedInterval);
		Utils.nonNull(readsSource);
		//Utils.validateArg(paddedInterval.contains(interval), "The padded interval must contain the un-padded interval");

		this.interval = interval;
		this.paddedInterval = paddedInterval;
		this.readsSource = readsSource;
	}

	/**
	 * Create a new Shard spanning the specified interval, with no additional
	 * padding
	 *
	 * @param interval
	 *            the genomic span covered by this shard
	 * @param readsSource
	 *            source of reads from which to populate this shard
	 */
	public LocalReadShard(final GenomeLocation interval, final ReadsDataSource readsSource) {
		this(interval, interval, readsSource);
	}

	/**
	 * Reads in this shard will be filtered using this filter before being
	 * returned. Read filtering will be performed before any requested
	 * downsampling.
	 *
	 * @param filter
	 *            filter to use (may be null, which signifies that no filtering
	 *            is to be performed)
	 */
	public void setReadFilter(final ReadFilter filter) {
		this.readFilter = filter;
	}

	/**
	 * Reads in this shard will be downsampled using this downsampler before
	 * being returned. Downsampling will be performed after any requested read
	 * filtering.
	 *
	 * @param downsampler
	 *            downsampler to use (may be null, which signifies that no
	 *            downsampling is to be performed)
	 */
	public void setDownsampler(final ReadsDownsampler downsampler) {
		this.downsampler = downsampler;
	}

	/**
	 * @return the interval this shard spans
	 */
	@Override
	public GenomeLocation getInterval() {
		return interval;
	}

	/**
	 * @return the interval this shard spans, potentially with additional
	 *         padding on each side
	 */
	@Override
	public GenomeLocation getPaddedInterval() {
		return paddedInterval;
	}

	/**
	 * @return number of bases of padding to the left of our interval
	 */
	public int numLeftPaddingBases() {
		return interval.getStart() - paddedInterval.getStart();
	}

	/**
	 * @return number of bases of padding to the right of our interval
	 */
	public int numRightPaddingBases() {
		return paddedInterval.getEnd() - interval.getEnd();
	}

	/**
	 * @param loc
	 *            Locatable to test
	 * @return true if loc is completely contained within this shard's interval,
	 *         otherwise false
	 */
	public boolean contains(final Locatable loc) {
		Utils.nonNull(loc);
		return interval.contains(loc);
	}

	/**
	 * @param loc
	 *            Locatable to test
	 * @return true if loc starts within this shard's interval, otherwise false
	 */
	public boolean containsStartPosition(final Locatable loc) {
		Utils.nonNull(loc);
		return interval.contains(new GenomeLocation(loc.getContig(), loc.getStart(), loc.getStart()));
	}

	/**
	 * @return an iterator over reads in this shard, as filtered using the
	 *         configured read filter and downsampled using the configured
	 *         downsampler; reads are lazily loaded rather than pre-loaded
	 *
	 *         Note that any read filtering is always performed before any
	 *         downsampling.
	 */
	@Override
	public Iterator<GaeaSamRecord> iterator() {
		readsSource.set(paddedInterval, readFilter, downsampler);
		
		return readsSource.iterator();
		//Iterator<GaeaSamRecord> readsIterator = readsSource.query(paddedInterval,readFilter);

		/*if (readFilter != null) {
			readsIterator = new ReadFilteringIterator(readsIterator, readFilter);
		}

		if (downsampler != null) {
			readsIterator = new ReadsDownsamplingIterator(readsIterator, downsampler);
		}

		return readsIterator;*/
	}

	/**
	 * @return a List containing all reads in this shard, pre-loaded, filtered
	 *         using the configured read filter, and downsampled using the
	 *         configured downsampler
	 *
	 *         Call {@link #iterator} instead to avoid pre-loading all reads at
	 *         once.
	 *
	 *         Note that any read filtering is always performed before any
	 *         downsampling.
	 */
	public List<GaeaSamRecord> loadAllReads() {
		return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED), false)
				.collect(Collectors.toList());
	}

	/**
	 * Divide an interval into LocalReadShards. Each shard will cover up to
	 * shardSize bases, include shardPadding bases of extra padding on either
	 * side, and begin shardSize bases after the previous shard (ie., shards
	 * will not overlap except potentially in the padded regions).
	 *
	 * @param interval
	 *            interval to shard; must be on the contig according to the
	 *            provided dictionary
	 * @param shardSize
	 *            desired shard size; intervals larger than this will be divided
	 *            into shards of up to this size
	 * @param shardPadding
	 *            desired shard padding; each shard's interval will be padded on
	 *            both sides by this number of bases (may be 0)
	 * @param readsSource
	 *            data source for reads
	 * @param dictionary
	 *            sequence dictionary for reads
	 * @return List of {@link LocalReadShard} objects spanning the interval
	 */
	public static List<LocalReadShard> divideIntervalIntoShards(final GenomeLocation interval, final int shardSize,
			final int shardPadding, final ReadsDataSource readsSource, final SAMSequenceDictionary dictionary) {
		return divideIntervalIntoShards(interval, shardSize, shardSize, shardPadding, readsSource, dictionary);
	}

	/**
	 * Divide an interval into LocalReadShards. Each shard will cover up to
	 * shardSize bases, include shardPadding bases of extra padding on either
	 * side, and begin shardStep bases after the previous shard.
	 *
	 * @param interval
	 *            interval to shard; must be on the contig according to the
	 *            provided dictionary
	 * @param shardSize
	 *            desired shard size; intervals larger than this will be divided
	 *            into shards of up to this size
	 * @param shardStep
	 *            each shard will begin this many bases away from the previous
	 *            shard
	 * @param shardPadding
	 *            desired shard padding; each shard's interval will be padded on
	 *            both sides by this number of bases (may be 0)
	 * @param readsSource
	 *            data source for reads
	 * @param dictionary
	 *            sequence dictionary for reads
	 * @return List of {@link LocalReadShard} objects spanning the interval
	 */
	public static List<LocalReadShard> divideIntervalIntoShards(final GenomeLocation interval, final int shardSize,
			final int shardStep, final int shardPadding, final ReadsDataSource readsSource,
			final SAMSequenceDictionary dictionary) {
		Utils.nonNull(readsSource);
		return Shard.divideIntervalIntoShards(interval, shardSize, shardStep, shardPadding, dictionary).stream()
				.map(shardBoundary -> new LocalReadShard(shardBoundary.getInterval(), shardBoundary.getPaddedInterval(),
						readsSource))
				.collect(Collectors.toList());
	}
}
