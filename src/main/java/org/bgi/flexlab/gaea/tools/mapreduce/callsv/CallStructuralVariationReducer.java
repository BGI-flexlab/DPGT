package org.bgi.flexlab.gaea.tools.mapreduce.callsv;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.NullWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;
import org.bgi.flexlab.gaea.tools.callsv.BuildConnection;
import org.bgi.flexlab.gaea.tools.callsv.LinkRegion;
import org.bgi.flexlab.gaea.tools.callsv.NewMapKey;
import org.bgi.flexlab.gaea.tools.callsv.Reads;
import org.bgi.flexlab.gaea.tools.callsv.Region;
import org.bgi.flexlab.gaea.tools.callsv.SamWritable;

public class CallStructuralVariationReducer extends Reducer<NewMapKey, SamWritable, NullWritable, Text>{

	/**
	 * 程序参数类Options，用于保存程序输入的参数
	 */
	private CallStructuralVariationOptions option = new CallStructuralVariationOptions();
	private Configuration conf;
	private BuildConnection bc;
	

	@Override
	protected void setup(Context context) throws IOException, InterruptedException {
		conf = context.getConfiguration();
		option.getOptionsFromHadoopConf(conf);
		bc = new BuildConnection(conf, option);
		bc.setUpperLower();
		
	}
	
	@Override
	protected void reduce(NewMapKey key, Iterable<SamWritable> values, Context context) throws IOException, InterruptedException {
		
		List<SamWritable> APRs = bc.getAPRs(key, values);
		if(APRs.isEmpty()) return;
			
		Map<Integer, Region> regInfoMap = bc.getRegion(APRs);
		if(regInfoMap.isEmpty()) return;
		
		Map<LinkRegion, List<Reads>> linkRegMap = bc.buildLink();
		if(linkRegMap.isEmpty()) return;
		
		for(Map.Entry<LinkRegion, List<Reads>> linkRegEntry : linkRegMap.entrySet()) {
			LinkRegion linkReg = linkRegEntry.getKey();
			List<Reads> reads = linkRegEntry.getValue();
			
			/**
			 *得到相互连通的两个区域，firstReg和secondReg
			 */
			Region firstReg = regInfoMap.get(linkReg.getFirstRegion());
			Region secondReg = regInfoMap.get(linkReg.getSecondRegion());
			
			Text sv = bc.svCaller(linkReg, reads, firstReg, secondReg);
			
			if(sv != null) {
				context.write(NullWritable.get(), sv);
			}
		}
		
		bc.close();
	}

}
