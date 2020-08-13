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
package org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin;

import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.model.GaussianMixtureModel;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.VariantDatum;
import org.bgi.flexlab.gaea.util.RandomUtils;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class VariantRecalibrationEngine {

	public final static double MIN_ACCEPTABLE_LOD_SCORE = -20000.0;

    // the unified argument collection
    final private VCFQualityControlOptions options;

    private final static double MIN_PROB_CONVERGENCE = 2E-2;

    /////////////////////////////
    // Public Methods to interface with the Engine
    /////////////////////////////

    public VariantRecalibrationEngine( final VCFQualityControlOptions options ) {
        this.options = options;
    }

    public GaussianMixtureModel generateModel( final List<VariantDatum> data ) {
        final GaussianMixtureModel model = new GaussianMixtureModel( options.getMaxGaussians(), data.get(0).annotations.length, options.getShrinkage(), options.getDirichletParamenter(), options.getPriorCounts() );
        variationalBayesExpectationMaximization( model, data );
        return model;
    }

    public void evaluateData( final List<VariantDatum> data, final GaussianMixtureModel model, final boolean evaluateContrastively ) {
        if( !model.isModelReadyForEvaluation ) {
            try {
                model.precomputeDenominatorForEvaluation();
            } catch( Exception e ) {
                model.failedToConverge = true;
                return;
            }
        }
        
        for( final VariantDatum datum : data ) {
            final double thisLod = evaluateDatum( datum, model );
            if( Double.isNaN(thisLod) ) {
                model.failedToConverge = true;
                return;
            }

            datum.lod = ( evaluateContrastively ?
                            ( Double.isInfinite(datum.lod) ? // positive model said negative infinity
                                    ( MIN_ACCEPTABLE_LOD_SCORE + RandomUtils.getRandomGenerator().nextDouble() * MIN_ACCEPTABLE_LOD_SCORE ) // Negative infinity lod values are possible when covariates are extremely far away from their tight Gaussians
                                    : datum.prior + datum.lod - thisLod) // contrastive evaluation: (prior + positive model - negative model)
                            : thisLod ); // positive model only so set the lod and return
        }
    }
    
    public void calculateWorstPerformingAnnotationMultiTask(final List<VariantDatum> data, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel,int taskNumber){
    	ExecutorService fixedThreadPool = Executors.newFixedThreadPool(taskNumber);
    	int size = data.size();
    	for(int i = 0 ; i < size ; i++){
    		final VariantDatum datum = data.get(i);
    		
    		fixedThreadPool.execute(new Runnable() {
				public void run() {
					int worstAnnotation = -1;
					double minProb = Double.MAX_VALUE;
					for( int iii = 0; iii < datum.annotations.length; iii++ ) {
						final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(datum, iii);
						final Double badProbLog10 = badModel.evaluateDatumInOneDimension(datum, iii);
						if( goodProbLog10 != null && badProbLog10 != null ) {
							final double prob = goodProbLog10 - badProbLog10;
							if(prob < minProb) { minProb = prob; worstAnnotation = iii; }
						}
					}
					datum.worstAnnotation = worstAnnotation;
			}});
        }
    }

    public void calculateWorstPerformingAnnotation( final List<VariantDatum> data, final GaussianMixtureModel goodModel, final GaussianMixtureModel badModel ) {
        for( final VariantDatum datum : data ) {
            int worstAnnotation = -1;
            double minProb = Double.MAX_VALUE;
            for( int iii = 0; iii < datum.annotations.length; iii++ ) {
                final Double goodProbLog10 = goodModel.evaluateDatumInOneDimension(datum, iii);
                final Double badProbLog10 = badModel.evaluateDatumInOneDimension(datum, iii);
                if( goodProbLog10 != null && badProbLog10 != null ) {
                    final double prob = goodProbLog10 - badProbLog10;
                    if(prob < minProb) { minProb = prob; worstAnnotation = iii; }
                }
            }
            datum.worstAnnotation = worstAnnotation;
        }
    }


    /////////////////////////////
    // Private Methods used for generating a GaussianMixtureModel
    /////////////////////////////

    private void variationalBayesExpectationMaximization( final GaussianMixtureModel model, final List<VariantDatum> data ) {
        model.initializeRandomModel( data, options.getNumKMeansIterations() );

        // The VBEM loop
        model.normalizePMixtureLog10();
        
        model.expectationStepMulti( data ,0);
    	
        double currentChangeInMixtureCoefficients;
        int iteration = 0;
        long start = System.currentTimeMillis();
        while( iteration < options.getMaxIterations() ) {
            iteration++;
            model.maximizationStepMulti( data);
            currentChangeInMixtureCoefficients = model.normalizePMixtureLog10();
            model.expectationStepMulti( data , iteration);
            if( iteration % 5 == 0 ) { // cut down on the number of output lines so that users can read the warning messages
            }
            if( iteration > 2 && currentChangeInMixtureCoefficients < MIN_PROB_CONVERGENCE ) {
                break;
            }
        }
        
        System.err.println("while time:"+(System.currentTimeMillis()-start)/1000+"s");
        
        model.evaluateFinalModelParameters( data );
    }

    /////////////////////////////
    // Private Methods used for evaluating data given a GaussianMixtureModel
    /////////////////////////////

    private double evaluateDatum( final VariantDatum datum, final GaussianMixtureModel model ) {
        return model.evaluateDatum( datum );
    }
}
