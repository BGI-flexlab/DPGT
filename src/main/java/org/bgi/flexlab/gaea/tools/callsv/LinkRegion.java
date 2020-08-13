package org.bgi.flexlab.gaea.tools.callsv;

import java.util.List;


public class LinkRegion implements Comparable<LinkRegion>{
	
	private int firstRegion;
	private int secondRegion;
	
	
	public LinkRegion(int firstRegion, int secondRegion) {
		this.firstRegion = Math.min(firstRegion, secondRegion);
		this.secondRegion = Math.max(firstRegion, secondRegion);
	}
	
	public LinkRegion(List<Integer> l) {
		this.firstRegion = Math.min(l.get(0), l.get(1));
		this.secondRegion = Math.max(l.get(0), l.get(1));
	}

	public int getFirstRegion() {
		return firstRegion;
	}

	public void setFirstRegion(int firstRegion) {
		this.firstRegion = firstRegion;
	}

	public int getSecondRegion() {
		return secondRegion;
	}

	public void setSecondRegion(int secondRegion) {
		this.secondRegion = secondRegion;
	}

	
	@Override
	public boolean equals(Object obj) {
		if(!(obj instanceof LinkRegion))
			throw new ClassCastException("can not cast to LinkRegion class");
		
		LinkRegion l = (LinkRegion)obj;
		return this.firstRegion==l.firstRegion && this.secondRegion==l.secondRegion;
	}

	@Override
	public String toString() {
		return firstRegion + "\t" + secondRegion;
	}

	public int compareTo(LinkRegion o) {
		if(this.firstRegion != o.firstRegion)
			return this.firstRegion - o.firstRegion;
		else
			return this.secondRegion - o.secondRegion;
		
	}


}
