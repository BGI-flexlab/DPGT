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
package org.bgi.flexlab.gaea.util;

import java.util.Random;

public class RandomUtils {
	public RandomUtils() {}

	private static final long RANDOM_SEED = 47382911L;
	private static Random randomGenerator = new Random(RANDOM_SEED);

	public static Random getRandomGenerator() {
		return randomGenerator;
	}

	public static Random getNewRandomGenerator() {
		return new Random(RANDOM_SEED);
	}

	public static void resetRandomGenerator() {
		randomGenerator.setSeed(RANDOM_SEED);
	}

	public static void resetRandomGenerator(long seed) {
		randomGenerator.setSeed(seed);
	}
}
