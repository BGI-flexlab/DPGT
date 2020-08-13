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
 * Copyright (C)  2016  Pablo Cingolani(pcingola@users.sourceforge.net)
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package org.bgi.flexlab.gaea.tools.annotator.util;

import java.util.Date;

public class Timer {

	private static Timer timer = new Timer(); // Keep track of time (since first class instantiation)

	boolean useMiliSecs = false;
	Date start;
	Date end;

	/**
	 * Show absolute timer value and a message
	 * @param msg
	 */
	public static void show(Object msg) {
		System.out.println(timer + "\t" + (msg == null ? "null" : msg.toString()));
	}

	/**
	 * Show absolute timer value and a message on STDERR
	 * @param msg
	 */
	public static void showStdErr(Object msg) {
		System.err.println(timer + "\t" + (msg == null ? "null" : msg.toString()));
	}

	public static String toString(long elapsedMs, boolean useMiliSecs) {
		long delta = elapsedMs;
		long days = delta / (24 * 60 * 60 * 1000);
		long hours = (delta % (24 * 60 * 60 * 1000)) / (60 * 60 * 1000);
		long mins = (delta % (60 * 60 * 1000)) / (60 * 1000);
		long secs = (delta % (60 * 1000)) / (1000);
		long ms = (delta % 1000);

		if (days > 0) {
			if (useMiliSecs) return String.format("%d days %02d:%02d:%02d.%03d", days, hours, mins, secs, ms);
			return String.format("%d days %02d:%02d:%02d", days, hours, mins, secs);
		}

		if (useMiliSecs) return String.format("%02d:%02d:%02d.%03d", hours, mins, secs, ms);
		return String.format("%02d:%02d:%02d", hours, mins, secs);
	}

	public Timer() {
		start = new Date();
	}

	/**
	 * Elapsed time in milliseconds
	 */
	public long elapsed() {
		if (end != null) return end.getTime() - start.getTime();
		Date now = new Date();
		return now.getTime() - start.getTime();
	}

	public void end() {
		end = new Date();
	}

	public void setUseMiliSecs(boolean useMiliSecs) {
		this.useMiliSecs = useMiliSecs;
	}

	public void start() {
		start = new Date();
		end = null;
	}

	@Override
	public String toString() {
		return toString(elapsed(), useMiliSecs);
	}

}
