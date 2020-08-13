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
 *
 * This file incorporates work covered by the following copyright and 
 * Permission notices:
 *
 * Copyright (c) 2009-2012 The Broad Institute
 *  
 *     Permission is hereby granted, free of charge, to any person
 *     obtaining a copy of this software and associated documentation
 *     files (the "Software"), to deal in the Software without
 *     restriction, including without limitation the rights to use,
 *     copy, modify, merge, publish, distribute, sublicense, and/or sell
 *     copies of the Software, and to permit persons to whom the
 *     Software is furnished to do so, subject to the following
 *     conditions:
 *  
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *  
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *     OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.realigner.event;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;

import java.util.ArrayList;

public class Event {
	public enum EVENT_TYPE {
		POINT_EVENT, INDEL_EVENT, BOTH
	}

	private int eventStartPosition = -1;
	private int eventStopPosition = -1;
	private int furthestStopPosition;

	private GenomeLocation location;
	private ArrayList<Integer> pointEvents = null;

	private EVENT_TYPE type;

	public Event(GenomeLocation location, int furthest, EVENT_TYPE type) {
		this.type = type;
		this.furthestStopPosition = furthest;
		this.location = location;

		if (type != EVENT_TYPE.POINT_EVENT) {
			this.eventStartPosition = location.getStart();
			this.eventStopPosition = location.getStop();
		}

		pointEvents = new ArrayList<Integer>();
		if (type != EVENT_TYPE.INDEL_EVENT) {
			pointEvents.add(location.getStart());
		}
	}
	
	public String toString(){
    	return location.getContig()+"-"+eventStartPosition+"-"+eventStopPosition+"->"+furthestStopPosition+":"+type.toString();
    }

	public GenomeLocation createLocation(GenomeLocationParser parser) {
		return parser.createGenomeLocation(location.getContig(),
				eventStartPosition, eventStopPosition);
	}

	public boolean isValidEvent(GenomeLocationParser parser, int maxIntervalSize) {
		return parser.isValidGenomeLocation(location.getContig(),
				eventStartPosition, eventStopPosition, true)
				&& (eventStopPosition - eventStartPosition < maxIntervalSize);
	}

	public GenomeLocation getLocation() {
		return this.location;
	}
	
	public GenomeLocation getLocation(GenomeLocationParser parser) {
		return parser.createGenomeLocation(location.getContig(), eventStartPosition, eventStopPosition);
	}

	public EVENT_TYPE getEventType() {
		return type;
	}

	public ArrayList<Integer> getPointEvents() {
		return this.pointEvents;
	}

	public int peek() {
		if (this.pointEvents.size() == 0)
			return -1;
		return this.pointEvents.get(0);
	}

	public int last() {
		if (this.pointEvents.size() == 0)
			return -1;
		return this.pointEvents.get(pointEvents.size() - 1);
	}

	public int getFurthestStop() {
		return this.furthestStopPosition;
	}

	public int getStop() {
		return this.eventStopPosition;
	}

	public int getStart() {
		return this.eventStartPosition;
	}

	public boolean canBeMerge(Event beMerge) {
		return getLocation().getContigIndex() == beMerge.getLocation()
				.getContigIndex()
				&& this.furthestStopPosition >= beMerge.getLocation()
						.getStart();
	}

	public void set(int start, int stop, int furthStop) {
		this.eventStartPosition = start;
		this.eventStopPosition = stop;
		this.furthestStopPosition = furthStop;
	}

	public void merge(Event event, int snpWin) {
		if (event.getEventType() != EVENT_TYPE.POINT_EVENT) {
			this.furthestStopPosition = event.getFurthestStop();
			this.eventStopPosition = event.getStop();
			if (this.eventStartPosition == -1)
				this.eventStartPosition = event.getStart();
		}
		if (event.getEventType() != EVENT_TYPE.INDEL_EVENT) {
			int newPosition = event.peek();
			
			if (pointEvents.size() > 0) {
				int lastPosition = last();
				if (newPosition - lastPosition < snpWin) {
					int minStart = eventStartPosition == -1 ? lastPosition
							: Math.min(eventStartPosition, lastPosition);
					set(minStart, Math.max(eventStopPosition, newPosition),
							event.getFurthestStop());
				} else if (eventStartPosition == -1
						&& event.getStart() != eventStartPosition) {
					set(event.getStart(), event.getStop(),
							event.getFurthestStop());
				}
			}
			pointEvents.add(newPosition);
		}
	}
}
