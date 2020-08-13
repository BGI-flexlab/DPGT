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
package org.bgi.flexlab.gaea.tools.mapreduce.realigner;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.io.Writable;
import org.bgi.flexlab.gaea.data.mapreduce.input.header.SamHdfsFileHeader;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.framework.tools.mapreduce.PairEndAggregatorMapper;
import org.bgi.flexlab.gaea.tools.recalibrator.report.RecalibratorReport;

public class FixmateMapper extends PairEndAggregatorMapper {
	private final String DefaultReadGroup = "UK";
	private String RG;
	private SAMFileHeader header = null;
	private Text readID = new Text();
	private RealignerExtendOptions option = new RealignerExtendOptions();
	private RecalibratorReport report = null;
	private SamRecordWritable writable = new SamRecordWritable();

	@Override
	public void setup(Context context) {
		Configuration conf = context.getConfiguration();
		header = SamHdfsFileHeader.getHeader(conf);
		option.getOptionsFromHadoopConf(conf);

		if (option.isRecalibration()) {
			String input = conf.get(Realigner.RECALIBRATOR_REPORT_TABLE_NAME);
			if (input == null)
				throw new RuntimeException("bqsr report table is null!!!");
			RecalibratorOptions bqsrOption = option.getBqsrOptions();
			report = new RecalibratorReport(input, header, 0,
					bqsrOption.PRESERVE_QSCORES_LESS_THAN);
		}
	}

	protected Writable getKey(Writable keyin, Writable valuein) {
		if(!option.isRealignment())
			return NullWritable.get();
		
		SAMRecord record = ((SamRecordWritable) valuein).get();
		GaeaSamRecord sam = new GaeaSamRecord(header, record);
		RG = (String) sam.getAttribute("RG");
		if (RG == null)
			RG = DefaultReadGroup;

		readID.set(RG + ":" + sam.getReadName());
		return readID;
	}

	protected Writable getValue(Writable value) {
		if (!option.isRecalibration())
			return value;

		if (value instanceof SamRecordWritable) {
			SamRecordWritable temp = (SamRecordWritable) value;

			GaeaSamRecord sam = new GaeaSamRecord(header, temp.get());
			report.readRecalibrator(sam);
			writable.set(sam);
			return writable;
		}

		return value;
	}
}
