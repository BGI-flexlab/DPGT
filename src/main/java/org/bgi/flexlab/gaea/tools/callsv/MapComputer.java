package org.bgi.flexlab.gaea.tools.callsv;

import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.bgi.flexlab.gaea.tools.mapreduce.callsv.CallStructuralVariationOptions;

import htsjdk.samtools.SAMRecord;

public class MapComputer {
	
	private Map<Integer, Integer> insertsize ;
	
	private CallStructuralVariationOptions option;
	
	public MapComputer() {
		insertsize = new TreeMap<Integer, Integer>();
		option = new CallStructuralVariationOptions();
	}

	public Map<Integer, Integer> getInsertsize() {
		return insertsize;
	}

	public void setInsertsize(Map<Integer, Integer> insertsize) {
		this.insertsize = insertsize;
	}
	
	public CallStructuralVariationOptions getOption() {
		return option;
	}

	public void setOption(CallStructuralVariationOptions option) {
		this.option = option;
	}

	/**
	 * saveInsert方法<br>
	 * 将正常的reads的insert size保存到中间文件中，中间文件名:/HDFSdir/Sort/LibConf/UUID <br>
	 * 一般一个map分片会选取100个正常的insert size<br>
	 * <br>
	 * @param record bam文件中每一个记录，也就是每一条read的比对情况
	 * @throws IOException 抛出IO异常
	 */
	public void saveInsert(SAMRecord record) {
		int insert = record.getInferredInsertSize();
		if(insert > 2000)
			return;
		if(insert <= 0)
			return;
		if(record.getMappingQuality() < option.getMinqual())
			return;
		if(!record.getReadPairedFlag())
			return;
		if(!record.getProperPairFlag())
			return;
		
		updateInsertsize(insert);
	}

	private void updateInsertsize(int insert) {
		int num = 0;
		if(!insertsize.containsKey(insert))
			num = 1;
		else {
			num = insertsize.get(insert);
			num ++ ;
		}
		insertsize.put(insert, num);
	}

	
	/**
	 * readClassify方法<br>
	 * 用于判断reads是否是异常reads（APRs），如果是，根据flag值判断属于哪一种类型；<br>
	 * 将reads信息保存到Format对象中，context.write到reducer中<br>
	 * @param record bam文件中每一个记录，也就是每一条read的比对情况
	 * @throws IOException 抛出IO异常
	 * @throws InterruptedException 抛出中断异常
	 */
	public MapContextWriter readClassify(SAMRecord record) {
		
		MapContextWriter res = new MapContextWriter();
		
		String type = null;
		int start1 = record.getAlignmentStart();
		int start2 = record.getMateAlignmentStart();
		
		if(record.getMappingQuality() <= option.getMinqual()) //low quality
			return null;
		else if(record.getDuplicateReadFlag())
		//else if(record.getDuplicateReadFlag() || record.getNotPrimaryAlignmentFlag()) //重复
			return null;
		else if(record.getReadPairedFlag()) { //pair read
			if(record.getReadUnmappedFlag()) //unmap
				return null;
			else if(record.getMateUnmappedFlag())  //mate unmap
				return null;
			else if(Math.abs(record.getInferredInsertSize()) >= option.getMaxsvsize()) //too long sv size
				return null;
			else if(!(record.getMateReferenceName().equals("=") || record.getMateReferenceName().equals(record.getReferenceName())))
				type = "Diff";
			else if(record.getReadNegativeStrandFlag() && record.getMateNegativeStrandFlag()) //--
				type = "FF_RR";
			else if(!record.getReadNegativeStrandFlag() && !record.getMateNegativeStrandFlag()) //++
				type = "FF_RR";
			else {
				
				if(start1 < start2) {
					if(record.getReadNegativeStrandFlag())
						type = "RF";
					else
						type = "FR";
				}else if(start1 >= start2){
					if(record.getReadNegativeStrandFlag())
						type = "FR";
					else
						type = "RF";
				}else
					return null;
			}
		}
		
		if(!type.equals(null) && !"Diff".equals(type)) {
			res.setAll(record, type);
		}else {
			return null;
		}
		return res;
	}
}
