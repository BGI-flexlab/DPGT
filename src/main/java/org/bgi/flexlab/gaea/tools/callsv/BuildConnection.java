package org.bgi.flexlab.gaea.tools.callsv;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.bgi.flexlab.gaea.tools.mapreduce.callsv.CallStructuralVariationOptions;


/**
 * BuildConnection类，用作CallSV主要的计算类
 * @author Huifang Lu
 *
 */
public class BuildConnection {
	
	private Configuration conf;
	
	/**
	 * 保存了程序的输入参数
	 */
	private CallStructuralVariationOptions options;
	
	private int mean;
	
	private float upper;
	
	private float lower;
	
	private float dist;
	
	private int ref_length;

	/**
	 * Map集合，key是ReadName，value是Reads的对象，保存了reads的信息
	 */
	private Map<String, Reads> readInfoMap;
	
	private BufferedReader br ;
	private FileSystem fs ;
	
	public BuildConnection() {
		this.conf = new Configuration();
		this.options = new CallStructuralVariationOptions();
		this.readInfoMap = new TreeMap<String, Reads>();
		this.mean = 0;
		this.upper = 0;
		this.lower = 0;
		this.dist = 0;
		this.ref_length = 0;
	}
	
	public BuildConnection(Configuration conf, CallStructuralVariationOptions options) {
		this.conf = conf;
		this.options = options;
		this.readInfoMap = new TreeMap<String, Reads>();
		this.mean = 0;
		this.upper = 0;
		this.lower = 0;
		this.dist = 0;
		this.ref_length = 0;
	}
	
	/**
	 * 带有参数para的构造函数<br>
	 * 如果使用此构造函数来构造对象，则初始化reg_id等变量
	 * @param para Parameter类型的参数，将上层的para传递进来
	 * @param dist Map集合，保存了每个染色体的dist和ref_length
	 */
	public BuildConnection(Configuration conf, CallStructuralVariationOptions options, int mean, float upper, float lower, float dist, int ref_length) {
		this.conf = conf;
		this.options = options;
		this.mean = mean;
		this.upper = upper;
		this.lower = lower;
		this.dist = dist;
		this.ref_length = ref_length;
		this.readInfoMap = new TreeMap<String, Reads>();
	}


	public Configuration getConf() {
		return conf;
	}


	public void setConf(Configuration conf) {
		this.conf = conf;
	}


	public CallStructuralVariationOptions getOptions() {
		return options;
	}


	public void setOptions(CallStructuralVariationOptions options) {
		this.options = options;
	}


	public float getDist() {
		return dist;
	}


	public void setDist(float dist) {
		this.dist = dist;
	}


	public int getRef_length() {
		return ref_length;
	}


	public void setRef_length(int ref_length) {
		this.ref_length = ref_length;
	}


	public Map<String, Reads> getReadInfoMap() {
		return readInfoMap;
	}

	public void setReadInfoMap(Map<String, Reads> readInfoMap) {
		this.readInfoMap = readInfoMap;
	}

	
	/**
	 * 将APRs划分区域，同时也保存了每一条reads的信息
	 * @param it 此参数输入的是一个iterator迭代器，保存了此reducer任务接收到的reads信息
	 * @return 一个Map集合，保存了所有reg信息，key是regid，value是Region对象
	 */
	public  Map<Integer, Region> getRegion(List<SamWritable> aprs) {
		Map<Integer, Region> regInfoMap = new TreeMap<Integer, Region>();
		
		Region reg = new Region();
		int regId = 0;
		
		for(SamWritable r : aprs) {
			/**
			 * chr不相同，并且间隔大于dist，满足划分为两个区域的条件，做break
			 * 但是还要判断是真的break还是要去掉的break
			 */
			if(!r.getChr().equals(reg.getChr()) || (r.getStart() - reg.getpStart()) > dist) { //break
				float coverage = reg.getRegCoverage();
				
				if(coverage > 0 && coverage < options.getMaxcoverage() &&  reg.getRegLength() > options.getMinlen()) { //real break
					regInfoMap.put(regId, reg);
					regId ++;
				}else { //false break
					for(String readid : reg.getRegReads()) {
						readInfoMap.remove(readid);  //delete false reads form readInfoMap
					}
					reg = null;
				}
				reg = new Region(r);
				reg.setRegId(regId);
			}
			
			/**
			 * 还是同一个区域，更新区域的信息
			 */
			reg.updateReg(r); // Update region
			regInfoMap.put(regId, reg);
			saveReadInfo(regId, r);
		}
		
		return regInfoMap;
	}

	/**
	 * 将read信息保存到readInfoMap集合中
	 * @param regId 当前区域编号
	 * @param r 当前read
	 */
	private void saveReadInfo(int regId, SamWritable r) {
		if (readInfoMap == null) {
			readInfoMap = new TreeMap<String, Reads>();
		}
		if(r == null) return;
		
		Reads read = readInfoMap.get(r.getReadName());
		if(read == null)
			read = new Reads(r);
		read.getReg().add(regId);
		readInfoMap.put(r.getReadName(), read);
		
	}
	
	

	/**
	 * 获取每一对有联系的区域，并保存其有联系的reads
	 * @return Map集合，key是有联系的两个区域，value是支持这两个区域有联系的reads列表
	 */
	public Map<LinkRegion, List<Reads>> buildLink() {
		Map<LinkRegion, List<Reads>> linkRegMap = new TreeMap<LinkRegion, List<Reads>>();
		
		for(Reads r: readInfoMap.values()) {
			/**
			 * 同一对reads不是比对到两个区域或者两个区域相等，去掉这对reads
			 */
			if(r.getReg().size() !=2)
				continue;
				
			/**
			 * 同一对reads只比对到两个区域，并且两个区域不相同，则这两个区域为相互连通的区域，保存下来
			 */
			LinkRegion tmpLinkReg = new LinkRegion(r.getReg());
			
			List<Reads> readList = linkRegMap.get(tmpLinkReg);
			if(readList == null)
				readList = new ArrayList<Reads>();
			readList.add(r);
			linkRegMap.put(tmpLinkReg, readList);
			
		}
		
		readInfoMap = null; 
		return linkRegMap;
	}
	
	
	public void setUpperLower() {
		Map<Integer, Integer> insert = readInsertFile(options.getHdfsdir() + "/Sort/LibConf/");
		
		int maxnum = 0;
		
		for(Map.Entry<Integer, Integer> entry : insert.entrySet()) {
			//get mean
			int num = entry.getValue();
			if(num > maxnum) {
				this.mean = entry.getKey();
				maxnum = num;
			}

		}

		long lowsum = 0;
		long upsum = 0;
		long lownum = 0;
		long upnum = 0;
		for(Map.Entry<Integer, Integer> entry : insert.entrySet()) {
			if(entry.getKey() < mean) {
				lowsum = (long) (lowsum + Math.pow((entry.getKey() - mean),2) * entry.getValue());
				lownum = lownum + entry.getValue();
			}else {
				upsum = (long) (upsum + Math.pow((entry.getKey() - mean),2) * entry.getValue());
				upnum = lownum + entry.getValue();
			}
		}
		
		float lowstd = (float) Math.sqrt(lowsum/(lownum-1));
		float upstd = (float) Math.sqrt(upsum/(upnum-1));
		
		upper = mean + options.getStdtimes() * upstd;
		lower = mean - options.getStdtimes() * lowstd;
				
	}
	
	
	public List<SamWritable> getAPRs(NewMapKey key, Iterable<SamWritable> values){
		
		List<SamWritable> aprs = new ArrayList<SamWritable>();
		float d = 100000000;
		int max_sd = 1000000000;
		int indel_num = 0;
		int min_pos = Integer.MAX_VALUE;
		int max_pos = 0;
		
		Iterator<SamWritable> vs = values.iterator();
		while (vs.hasNext()) {
			SamWritable value = vs.next();
			
			float tmp1 = mean - value.getReadLen()*2;
			d = Math.min(d, tmp1);
			d = Math.max(d, 50);
			
			min_pos = Math.min(value.getStart(), min_pos);
			max_pos = Math.max(value.getEnd(), max_pos);
			
			if(value.getInsert() > max_sd) {continue;}
			
			/**
			 * 判断FR类型的reads是DEL或者INS
			 */
			if(value.getType().equals("FR")) {
				if(Math.abs(value.getInsert()) > upper) {
					value.setType("FR_long");
					indel_num++;
				}else if(Math.abs(value.getInsert()) < lower) {
					value.setType("FR_short");
					indel_num++;
				}else {
					continue;
				}
			}
			
			value.setType(changeType(value.getType()));
			SamWritable f = new SamWritable(value.toString());
			aprs.add(f);
		}
		
		ref_length = max_pos - min_pos + 1;
		if(indel_num == 0 || ref_length == 0) {
			dist = d;
		}else {
			dist = Math.min(d, ref_length/indel_num);
		}
		dist = dist < 50 ? 50 : dist;

		return aprs;
	}

	
	/**
	 * 将reads的APRs所支持的类型转变成SV类型
	 * @param t reads的类型（"FR_long"， "FR_short"， "FF_RR"， "RF"， "Diff"）
	 * @return 返回字符串类型的SV Type（"DEL"， "INS"， "INV"， "ITX"， "CTX"）
	 */
	private String changeType(String t) {
		if(t.equals("FR_long"))
			return "DEL";
		else if(t.equals("FR_short"))
			return "INS";
		else if(t.equals("RF"))
			return "ITX";
		else if(t.equals("FF_RR"))
			return "INV";
		else if(t.equals("Diff"))
			return "CTX";
		else
			return null;

	}
	

	
	/**
	 * 遍历每一对有联系的区域，做最终的calling
	 * @param LinkRegion linkReg, List<Reads> reads, Region firstReg, Region secondReg
	 * @throws InterruptedException 抛出中断异常
	 * @throws IOException 抛出IO异常
	 */
	public Text svCaller(LinkRegion linkReg, List<Reads> reads, Region firstReg, Region secondReg) {
		
		/**
		 * 遍历这对相互连通对区域中的reads，保存每一种sv类型的信息
		 */
		Map<String, LinkRegType> linkRegTypeMap = saveTypeInfo(reads);
		
		/**
		 * 选取最终的sv类型finalType
		 */
		String finalType = selectFinalType(linkRegTypeMap);
		if(finalType == null)
			return null;
		
		/**
		 * 当最终的sv类型不为null时，计算finalType的分数score
		 */
		LinkRegType finalTypeInfo = linkRegTypeMap.get(finalType);
		int score = computeProbScore(firstReg, secondReg, finalTypeInfo);
		
		/**
		 * 当分数大于输出的分数时，输出到parts
		 */
		int size = finalTypeInfo.getSize()/finalTypeInfo.getReadNum();
		
		String writer = firstReg.firstToString() + "\t" + 
				secondReg.secondToString() + "\t" +
				finalType + "\t" + size + "\t" + 
				score + "\t" + finalTypeInfo.getReadNum();
		
		return new Text(writer);
	}
	
	
	
	/**
	 * 计算选定的最终类型的分数
	 * @param firstReg 这个SV相连的两个区域中第一个区域
	 * @param secondReg 这个SV相连的两个区域中第二个区域
	 * @param finalReg LinkRegType对象，保存最终类型的信息
	 * @return 计算得到的分数
	 */
	private int  computeProbScore(Region firstReg, Region secondReg, LinkRegType finalTypeInfo) {
		
		int totalRegSize = firstReg.getRegLength() + secondReg.getRegLength();
		
		double logP = 0;
		for(Integer libNum : finalTypeInfo.getLibNum().values()) {
			
			double lambda = (double)totalRegSize*libNum/ref_length; //total_reg_size*lib_read_num/ref_length
			Score sc = new Score();
			logP = logP + sc.logPoissionTailProb(libNum, lambda);
		}
		
		double phredQ = -10*(logP/Math.log(10));
		int score = (int) ((phredQ > 99) ? 99 : phredQ + 0.5);
		return score;
	}

	
	/**
	 * 遍历一对连通的区域中所有的reads，保存每一种Type的信息
	 * @param linkRegReads 一对连通的区域中所有的reads列表
	 * @return Map集合，key是type，value是LinkRegType对象
	 */
	private Map<String, LinkRegType> saveTypeInfo( List<Reads> linkRegReads) {
		
		Map<String, LinkRegType> linkRegType = new TreeMap<String, LinkRegType>();
		
		for(Reads r : linkRegReads) {
			LinkRegType typeInfo = linkRegType.get(r.getType());
			if(typeInfo==null)
				typeInfo = new LinkRegType(r);
			
			int size = Math.abs(r.getInsert()-mean);
			
			typeInfo.updateType(r, size);
			linkRegType.put(r.getType(), typeInfo);
		}
		return linkRegType;
	}

	
	/**
	 * 根据保存的每一种Type的信息，选取终的type
	 * @param linkRegType Map集合，保存了同一对连通的区域中，每一种类型的信息
	 * @return 最终选定的类型
	 */
	private String selectFinalType(Map<String, LinkRegType> linkRegType) {
		int finalNum = 0;
		String finalType = null;
		for(Map.Entry<String, LinkRegType> linkTypeEntry: linkRegType.entrySet()) {
			if(finalNum < linkTypeEntry.getValue().getReadNum()) {
				finalNum = linkTypeEntry.getValue().getReadNum();
				finalType = linkTypeEntry.getValue().getType();
			}
		}
		finalType = (finalNum >= options.getMinpair()) ? finalType : null;
		return finalType;
	}
	
	
	private Map<Integer, Integer> readInsertFile(String libconftxt){
		br = null; fs = null;
		
		Map<Integer, Integer> map = new TreeMap<Integer, Integer>();
		
		try {
			fs = FileSystem.get(this.conf);
			FileStatus[] flist = fs.listStatus(new Path(libconftxt));
			
			for(FileStatus file : flist) {
				
				FSDataInputStream fsopen = fs.open(file.getPath());
				br = new BufferedReader(new InputStreamReader(fsopen));
				
				String line = null;
				while((line = br.readLine())!= null) {
					String[] lines = line.split("\\t");
					
					int num;
					if(!map.containsKey(Integer.parseInt(lines[0]))) {
						num = Integer.parseInt(lines[1]);
					}else {
						int pnum = map.get(Integer.parseInt(lines[0]));
						num = pnum + Integer.parseInt(lines[1]);
					}
					map.put(Integer.parseInt(lines[0]), num);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return map;
	}
	
	
	public void close() {
		
		if(fs != null) {
			try {
				fs.close();
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		if(br != null) {
			try {
				br.close();
			}catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}
}
