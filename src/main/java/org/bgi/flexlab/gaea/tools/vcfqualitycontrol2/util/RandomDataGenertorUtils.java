package org.bgi.flexlab.gaea.tools.vcfqualitycontrol2.util;

import java.util.Random;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;

public class RandomDataGenertorUtils {
	 private static final long GATK_RANDOM_SEED = 47382911L;
	    private static final Random randomGenerator = new Random(GATK_RANDOM_SEED);
	    private static final RandomDataGenerator randomDataGenerator = new RandomDataGenerator(new Well19937c(GATK_RANDOM_SEED));

	    public static Random getRandomGenerator() { return randomGenerator; }
	    public static RandomDataGenerator getRandomDataGenerator() { return randomDataGenerator; }

	    public static void resetRandomGenerator() {
	        randomGenerator.setSeed(GATK_RANDOM_SEED);
	        randomDataGenerator.reSeed(GATK_RANDOM_SEED);
	    }
}
