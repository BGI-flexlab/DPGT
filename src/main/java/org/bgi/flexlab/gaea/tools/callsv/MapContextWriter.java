package org.bgi.flexlab.gaea.tools.callsv;

import htsjdk.samtools.SAMRecord;

public class MapContextWriter {
	
	private NewMapKey key;
	
	private SamWritable sam;
	
	public MapContextWriter() {
		key = new NewMapKey();
		sam = new SamWritable();
	}
	
	public NewMapKey getKey() {
		return key;
	}
	
	public void setKey(NewMapKey key) {
		this.key = key;
	}
	
	public SamWritable getSam() {
		return sam;
	}
	
	public void setSam(SamWritable sam) {
		this.sam = sam;
	}

	public void setAll (SAMRecord record, String type) {
		sam.set(record, type);
		key.setChr(record.getReferenceName());
		key.setPos(record.getAlignmentStart());
		key.setEnd(record.getAlignmentEnd());
	}
	
}
