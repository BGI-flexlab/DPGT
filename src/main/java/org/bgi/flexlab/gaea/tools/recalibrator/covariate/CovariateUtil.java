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
package org.bgi.flexlab.gaea.tools.recalibrator.covariate;

import htsjdk.samtools.SAMFileHeader;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.util.Pair;

import java.util.ArrayList;
import java.util.List;

public class CovariateUtil {
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private static List<Class<? extends Covariate>> initConvariate() {
		List<Class<? extends Covariate>> covariates = new ArrayList();
		covariates.add(ContextCovariate.class);
		covariates.add(CycleCovariate.class);
		covariates.add(QualityCovariate.class);
		covariates.add(ReadGroupCovariate.class);
		return covariates;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	private static List<Class<? extends RequiredCovariate>> initRequiredCovariate() {
		List<Class<? extends RequiredCovariate>> require = new ArrayList();
		require.add(QualityCovariate.class);
		require.add(ReadGroupCovariate.class);
		return require;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	private static List<Class<? extends OptionalCovariate>> initOptionalCovariate(boolean cycle) {
		List<Class<? extends OptionalCovariate>> optional = new ArrayList();
		optional.add(ContextCovariate.class);
		if(cycle)
			optional.add(CycleCovariate.class);
		return optional;
	}

	private static ArrayList<Covariate> addRequiredCovariatesToList(List<Class<? extends RequiredCovariate>> classes) {
		ArrayList<Covariate> dest = new ArrayList<Covariate>(classes.size());
		if (classes.size() != 2)
			throw new UserException("Require covariate had changed?");

		dest.add(new ReadGroupCovariate());
		dest.add(new QualityCovariate());
		return dest;
	}

	private static ArrayList<Covariate> addOptionalCovariatesToList(List<Class<? extends OptionalCovariate>> classes) {
		ArrayList<Covariate> dest = new ArrayList<Covariate>(classes.size());
		for (Class<?> covClass : classes) {
			try {
				final Covariate covariate = (Covariate) covClass.newInstance();
				dest.add(covariate);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return dest;
	}

	private static Pair<ArrayList<Covariate>, ArrayList<Covariate>> initializeCovariates(
			RecalibratorOptions option) {
		final List<Class<? extends Covariate>> covariateClasses = initConvariate();
		final List<Class<? extends RequiredCovariate>> requiredClasses = initRequiredCovariate();
		final List<Class<? extends OptionalCovariate>> standardClasses = initOptionalCovariate(false);

		final ArrayList<Covariate> requiredCovariates = addRequiredCovariatesToList(requiredClasses);
		ArrayList<Covariate> optionalCovariates = new ArrayList<Covariate>();
		if (!option.DO_NOT_USE_STANDARD_COVARIATES)
			optionalCovariates = addOptionalCovariatesToList(standardClasses);

		if (option.COVARIATES != null) {
			for (String requestedCovariateString : option.COVARIATES) {
				boolean foundClass = false;
				for (Class<? extends Covariate> covClass : covariateClasses) {
					if (requestedCovariateString.equalsIgnoreCase(covClass.getSimpleName())) {
						foundClass = true;
						if (!requiredClasses.contains(covClass)
								&& (option.DO_NOT_USE_STANDARD_COVARIATES || !standardClasses.contains(covClass))) {
							try {
								final Covariate covariate = covClass.newInstance();
								optionalCovariates.add(covariate);
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}
				}

				if (!foundClass) {
					throw new UserException.CommandLineException("The requested covariate type ("
							+ requestedCovariateString + ") isn't a valid covariate option.");
				}
			}
		}
		return new Pair<ArrayList<Covariate>, ArrayList<Covariate>>(requiredCovariates, optionalCovariates);
	}

	public static Covariate[] initializeCovariates(RecalibratorOptions option, SAMFileHeader mHeader) {
		Pair<ArrayList<Covariate>, ArrayList<Covariate>> covariates = initializeCovariates(option);
		ArrayList<Covariate> requiredCovariates = covariates.getFirst();
		ArrayList<Covariate> optionalCovariates = covariates.getSecond();

		Covariate[] requestedCovariates = new Covariate[requiredCovariates.size() + optionalCovariates.size()];
		int covariateIndex = 0;
		for (final Covariate covariate : requiredCovariates)
			requestedCovariates[covariateIndex++] = covariate;
		for (final Covariate covariate : optionalCovariates)
			requestedCovariates[covariateIndex++] = covariate;

		for (Covariate cov : requestedCovariates) {
			cov.initialize(option);
			if (cov instanceof ReadGroupCovariate)
				((ReadGroupCovariate) cov).initializeReadGroup(mHeader);
		}

		return requestedCovariates;
	}
}
