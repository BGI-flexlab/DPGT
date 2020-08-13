package org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata;

import java.util.HashMap;

import org.bgi.flexlab.gaea.data.exception.UserException;

public class ResourceTags {
	private HashMap<String,String> tags = null;
	
	public ResourceTags(){
		tags = new HashMap<String,String>();
		
		for(ResourceTag tag : ResourceTag.values()){
			tags.put(tag.tag, tag.property);
		}
	}
	
	public void parseTag(String resource) {
		String[] array = resource.split("#");
		for(String tag: array) {
			String[] keyValue = tag.split("=");
			if(!tags.containsKey(keyValue[0].toLowerCase()))
				throw new UserException("resource tags not contains key "+keyValue[0].toLowerCase());
			tags.put(keyValue[0].toLowerCase(), keyValue[1]);
		}
	}
	
	public String getProperty(ResourceTag tag) {
		return tags.get(tag.tag);
	}
}
