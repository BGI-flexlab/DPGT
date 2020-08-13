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
package org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter;

import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.Interval;
import org.bgi.flexlab.gaea.tools.bamqualtiycontrol.counter.CounterProperty.ReadType;

public class ReadsCounter {
	
	private long count;
	
	private ReadType type;
	
	private Interval interval;
	
	private CounterProperty[] properties;

	public ReadsCounter(CounterProperty... counterProperties) {
		properties = counterProperties;
		if(properties.length == 1) {
			type = (ReadType) properties[0];
		} else {
			type = (ReadType) properties[0];
			interval = (Interval) properties[1];
		}
	}
	
	public void update(ReadType type) {
		if(this.type == type)
			if(this.type != ReadType.PE)
				count++;
			else 
				count+=2;
	}
	
	public void update(ReadType type, Interval interval) {
		if(this.type == type && this.interval == interval)
			count++;
	}
	
	public String formatKey() {
		String key = "";
		for(CounterProperty property : properties)
			key += property.toString();
		return key;
	}
	
	public void setReadsCount(long count) {
		this.count += count;
	}
	
	public long getReadsCount() {
		return count;
	}
}
