package org.bgi.flexlab.gaea.tools.mapreduce.callsv;

import java.io.IOException;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataOutputStream;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.Mapper;
import org.bgi.flexlab.gaea.data.mapreduce.writable.SamRecordWritable;
import org.bgi.flexlab.gaea.tools.callsv.MapComputer;
import org.bgi.flexlab.gaea.tools.callsv.MapContextWriter;
import org.bgi.flexlab.gaea.tools.callsv.NewMapKey;
import org.bgi.flexlab.gaea.tools.callsv.SamWritable;
import org.seqdoop.hadoop_bam.FileVirtualSplit;

import htsjdk.samtools.SAMRecord;

public class CallStructuralVariationMapper extends Mapper<LongWritable, SamRecordWritable, NewMapKey, SamWritable>{

	private Configuration conf;
	private FSDataOutputStream out;
	private CallStructuralVariationOptions option = new CallStructuralVariationOptions();
	private MapComputer mc;
	

	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		conf = context.getConfiguration();
		option.getOptionsFromHadoopConf(conf);
		
		FileVirtualSplit input = (FileVirtualSplit)context.getInputSplit();
		String filename = input.getPath().getName();
		String libpath = option.getHdfsdir() + "/Sort/LibConf/" + filename + "-" + input.getStartVirtualOffset();
		out = FileSystem.get(conf).create(new Path(libpath));
		
		mc = new MapComputer();
		mc.setOption(option);
	}
	
	@Override
	protected void map(LongWritable key, SamRecordWritable value, Context context) throws IOException, InterruptedException {
		SAMRecord record = value.get();
		
		mc.saveInsert(record); //save insert
		MapContextWriter res = mc.readClassify(record);	//classify all reads
		
		if (res != null) {
			context.write(res.getKey(), res.getSam());
		}
		
	}

	@Override
	protected void cleanup(Context context) throws IOException, InterruptedException {
		
		Map<Integer, Integer> insertsize = mc.getInsertsize();
		for(Map.Entry<Integer, Integer> entry : insertsize.entrySet()) {
			String writer = entry.getKey() + "\t" + entry.getValue() + "\n";
			out.write(writer.getBytes());
		}
		out.flush();
		out.close();
		conf = null;
		option = null;
	}

}
