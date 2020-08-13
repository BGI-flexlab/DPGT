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
package org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.tranche;

import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.VariantDatum;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TrancheManager {
	 public static abstract class SelectionMetric {
	        String name = null;

	        public SelectionMetric(String name) {
	            this.name = name;
	        }

	        public String getName() { return name; }

	        public abstract double getThreshold(double tranche);
	        public abstract double getTarget();
	        public abstract void calculateRunningMetric(List<VariantDatum> data);
	        public abstract double getRunningMetric(int i);
	        public abstract int datumValue(VariantDatum d);
	    }

	    public static class NovelTiTvMetric extends SelectionMetric {
	        double[] runningTiTv;
	        double targetTiTv = 0;

	        public NovelTiTvMetric(double target) {
	            super("NovelTiTv");
	            targetTiTv = target; // compute the desired TiTv
	        }

	        public double getThreshold(double tranche) {
	            return fdrToTiTv(tranche, targetTiTv);
	        }

	        public double getTarget() { return targetTiTv; }

	        public void calculateRunningMetric(List<VariantDatum> data) {
	            int ti = 0, tv = 0;
	            runningTiTv = new double[data.size()];

	            for ( int i = data.size() - 1; i >= 0; i-- ) {
	                VariantDatum datum = data.get(i);
	                if ( ! datum.isKnown ) {
	                    if ( datum.isTransition ) { ti++; } else { tv++; }
	                    runningTiTv[i] = ti / Math.max(1.0 * tv, 1.0);
	                }
	            }
	        }

	        public double getRunningMetric(int i) {
	            return runningTiTv[i];
	        }

	        public int datumValue(VariantDatum d) {
	            return d.isTransition ? 1 : 0;
	        }
	    }

	    public static class TruthSensitivityMetric extends SelectionMetric {
	        double[] runningSensitivity;
	        int nTrueSites = 0;

	        public TruthSensitivityMetric(int nTrueSites) {
	            super("TruthSensitivity");
	            this.nTrueSites = nTrueSites;
	        }

	        public double getThreshold(double tranche) {
	            return 1.0 - tranche/100.0; // tranche of 1 => 99% sensitivity target
	        }

	        public double getTarget() { return 1.0; }

	        public void calculateRunningMetric(List<VariantDatum> data) {
	            int nCalledAtTruth = 0;
	            runningSensitivity = new double[data.size()];

	            for ( int i = data.size() - 1; i >= 0; i-- ) {
	                VariantDatum datum = data.get(i);
	                nCalledAtTruth += datum.atTruthSite ? 1 : 0;
	                runningSensitivity[i] = 1 - nCalledAtTruth / (1.0 * nTrueSites);
	            }
	        }

	        public double getRunningMetric(int i) {
	            return runningSensitivity[i];
	        }

	        public int datumValue(VariantDatum d) {
	            return d.atTruthSite ? 1 : 0;
	        }
	    }

	    public static List<Tranche> findTranches( final ArrayList<VariantDatum> data, final double[] tranches, final SelectionMetric metric, final VCFQualityControlOptions.Mode model ) {
	        return findTranches( data, tranches, metric, model, null );
	    }

	    public static List<Tranche> findTranches( final ArrayList<VariantDatum> data, final double[] trancheThresholds, final SelectionMetric metric, final VCFQualityControlOptions.Mode model, final File debugFile ) {

	        Collections.sort( data, new VariantDatum.VariantDatumLODComparator() );
	        metric.calculateRunningMetric(data);

	        if ( debugFile != null) { writeTranchesDebuggingInfo(debugFile, data, metric); }

	        List<Tranche> tranches = new ArrayList<Tranche>();
	        for ( double trancheThreshold : trancheThresholds ) {
	            Tranche t = findTranche(data, metric, trancheThreshold, model, tranches);

	            if ( t == null ) {
	                if ( tranches.size() == 0 )
	                    throw new UserException(String.format("Couldn't find any tranche containing variants with a %s > %.2f. Are you sure the truth files contain unfiltered variants which overlap the input data?", metric.getName(), metric.getThreshold(trancheThreshold)));
	                break;
	            }

	            tranches.add(t);
	        }

	        return tranches;
	    }

	    private static void writeTranchesDebuggingInfo(File f, List<VariantDatum> tranchesData, SelectionMetric metric ) {
	        try {
	            PrintStream out = new PrintStream(f);
	            out.println("Qual metricValue runningValue");
	            for ( int i = 0; i < tranchesData.size(); i++ ) {
	                VariantDatum  d = tranchesData.get(i);
	                int score = metric.datumValue(d);
	                double runningValue = metric.getRunningMetric(i);
	                out.printf("%.4f %d %.4f%n", d.lod, score, runningValue);
	            }
	            out.close();
	        } catch (FileNotFoundException e) {
	            throw new UserException.CouldNotCreateOutputFile(f, e);
	        }
	    }

	    public static Tranche findTranche( final List<VariantDatum> data, final SelectionMetric metric, final double trancheThreshold, final VCFQualityControlOptions.Mode model, List<Tranche> tranches ) {

	        double metricThreshold = metric.getThreshold(trancheThreshold);
	        int n = data.size();
	        for (int i = 0 ; i < n; i++ ) {
	            if ( metric.getRunningMetric(i) >= metricThreshold ) {
	                // we've found the largest group of variants with sensitivity >= our target truth sensitivity
	                Tranche t = trancheOfVariants(data, i, trancheThreshold, model, tranches);
	                return t;
	            }
	        }

	        return null;
	    }

	    public static Tranche trancheOfVariants( final List<VariantDatum> data, int minI, double ts, final VCFQualityControlOptions.Mode model, List<Tranche> tranches ) {
	        Tranche lastTranche = tranches.size() == 0 ? null : tranches.get(tranches.size() - 1);
	        int knownTi = 0, knownTv = 0, novelTi = 0, novelTv = 0;
	        if(lastTranche != null) {
	        	knownTi = lastTranche.knownTi;
	        	knownTv = lastTranche.knownTv;
	        	novelTi = lastTranche.novelTi;
	        	novelTv = lastTranche.novelTv;
	        }
	       
	        double minLod = data.get(minI).lod;
	        for ( final VariantDatum datum : data ) {
	            if ( datum.lod >= minLod ) {
	                //if( ! datum.isKnown ) System.out.println(datum.pos);
	                if ( datum.isKnown ) {
	                    if( datum.isSNP ) {
	                        if ( datum.isTransition ) { knownTi++; } else { knownTv++; }
	                    }
	                } else {
	                    if( datum.isSNP ) {
	                        if ( datum.isTransition ) { novelTi++; } else { novelTv++; }
	                    }
	                }
	            }
	            //data.remove(datum);
	        }

	        double knownTiTv = knownTi / Math.max(1.0 * knownTv, 1.0);
	        double novelTiTv = novelTi / Math.max(1.0 * novelTv, 1.0);

	        int accessibleTruthSites = countCallsAtTruth(data, Double.NEGATIVE_INFINITY);
	        int nCallsAtTruth = countCallsAtTruth(data, minLod);

	        return new Tranche(ts, minLod, knownTi, knownTv, knownTiTv, novelTi, novelTv, novelTiTv, accessibleTruthSites, nCallsAtTruth, model);
	    }

	    public static double fdrToTiTv(double desiredFDR, double targetTiTv) {
	            return (1.0 - desiredFDR / 100.0) * (targetTiTv - 0.5) + 0.5;
	    }

	    public static int countCallsAtTruth(final List<VariantDatum> data, double minLOD ) {
	        int n = 0;
	        for ( VariantDatum d : data ) { n += (d.atTruthSite && d.lod >= minLOD ? 1 : 0); }
	        return n;
	    }
}
