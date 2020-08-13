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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.bgi.flexlab.gaea.data.exception.UserException;
import org.bgi.flexlab.gaea.data.structure.header.VCFConstants;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControl;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.model.GaussianMixtureModel;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.model.VariantDataManager;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.VariantDatum;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata.VariantDatumMessenger;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.tranche.Tranche;
import org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.tranche.TrancheManager;
import org.bgi.flexlab.gaea.util.ExpandingArrayList;

import java.io.IOException;
import java.util.*;

public class VCFRecalibrationTable {
	/**
	 * parameter
	 */
    private VCFQualityControlOptions options; 
    
    /**
     * engine that generate recal table
     */
    private VariantRecalibrationEngine engine;
    
    /**
     * data Manager
     */
    private VariantDataManager dataManager;
    
    /**
     * tranches of recal data
     */
    private List<Tranche> tranches;
    
    
    /**
     * create dummy alleles to be used
     */
    private final List<Allele> alleles = new ArrayList<Allele>(2);   
    
    /**
     * to be used for the important INFO tags
     */
    private final HashMap<String, Object> attributes = new HashMap<String, Object>(3);
    
    /**
     * data of data manager
     */
    ExpandingArrayList<VariantDatum> data;
	Map<String, Map<Integer, ArrayList<Integer>>> dataIndex;

	/**
	 * construct function 
	 * @param options
	 * @throws IOException
	 */
	public VCFRecalibrationTable(VCFQualityControlOptions options) throws IOException {
		this.options = options;
    	engine = new VariantRecalibrationEngine(options);
		dataManager = new VariantDataManager( options.getUseAnnotations(), options );
		data = dataManager.getData();
		alleles.add(Allele.create("N", true));
	    alleles.add(Allele.create("<VQSR>", false)); 
	}
	
	/**
	 * add data
	 * @param data
	 */
	public void addData(VariantDatumMessenger data) {
		dataManager.addData(data);
	}
	
	/**
	 * generate recal table
	 * @return
	 */
	public final void getRecalibrationTable() {
		dataManager.normalizeData(); // Each data point is now (x - mean) / standard deviation
	
		long start = System.currentTimeMillis();
		
	    // Generate the positive model using the training data and evaluate each variant
	    final GaussianMixtureModel goodModel = engine.generateModel( dataManager.getTrainingData() );
	    engine.evaluateData( dataManager.getData(), goodModel, false );
	    
	    // Generate the negative model using the worst performing data and evaluate each variant contrastively
	  	final ExpandingArrayList<VariantDatum> negativeTrainingData = dataManager.selectWorstVariants( options.getPercentBadVariants(), options.getMinNumBadVariants() );
	 	GaussianMixtureModel badModel = engine.generateModel( negativeTrainingData );
	  	engine.evaluateData( dataManager.getData(), badModel, true );

    	// Detect if the negative model failed to converge because of too few points and/or too many Gaussians and try again
	  	while( badModel.failedToConverge && options.getMaxGaussians() > 4 ) {
	  		System.out.println("Negative model failed to converge. Retrying...");
	       	options.setMaxGaussians(options.getMaxGaussians() - 1);
	       	badModel = engine.generateModel( negativeTrainingData );
	      	engine.evaluateData( dataManager.getData(), goodModel, false );
	       	engine.evaluateData( dataManager.getData(), badModel, true );
	   	}
	  	
	  	System.err.println("generateModel time:"+(System.currentTimeMillis()-start)/1000+"s");
	
	   	if( badModel.failedToConverge || goodModel.failedToConverge ) {
	       	throw new UserException("NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example)");
	   	}

	   	//engine.calculateWorstPerformingAnnotation( dataManager.getData(), goodModel, badModel );
	   	engine.calculateWorstPerformingAnnotationMultiTask(negativeTrainingData, goodModel, badModel, options.getMaxGaussians());
	
	   	// Find the VQSLOD cutoff values which correspond to the various tranches of calls requested by the user
	   	final int nCallsAtTruth = TrancheManager.countCallsAtTruth( dataManager.getData(), Double.NEGATIVE_INFINITY );
	   	final TrancheManager.SelectionMetric metric = new TrancheManager.TruthSensitivityMetric( nCallsAtTruth );

	   	tranches = TrancheManager.findTranches( dataManager.getData(), options.getTsTranchesDouble(), metric, options.getMode() );
    	
	   	Collections.sort( tranches, new Tranche.TrancheTruthSensitivityComparator() );
	   	Tranche prev = null;
        for ( Tranche t : tranches ) {
            t.name = String.format("VQSRTranche%s%.2fto%.2f", t.model.toString(),(prev == null ? 0.0 : prev.ts), t.ts, t.model.toString());
            prev = t;
        }
	}
    
	public VariantContext getData(String chrName, int start, int end) {
    	ArrayList<Integer> index = dataIndex.get(chrName).get(start);
    	if(index == null)
    		return null;
    	for(int i : index) {
    		VariantDatum datum = data.get(i);
    		if(datum.loc.getStart() == start && datum.loc.getStop() == end) {
    			attributes.put(VCFConstants.END_KEY, datum.loc.getStop());
                attributes.put(VCFQualityControl.VQS_LOD_KEY, String.format("%.4f", datum.lod));
                attributes.put(VCFQualityControl.CULPRIT_KEY, (datum.worstAnnotation != -1 ? options.getUseAnnotations().get(datum.worstAnnotation) : "NULL"));

                VariantContextBuilder builder = new VariantContextBuilder("VQSR", datum.loc.getContig(), datum.loc.getStart(), datum.loc.getStop(), alleles).attributes(attributes);
                return builder.make();
    		}
    		else if(datum.loc.getStart() > start)
    			break;
    	}
    	return null;
    }
		
    public void indexData() {
    	dataIndex = new HashMap<String, Map<Integer,ArrayList<Integer>>>();
    	int i = 0;
    	while(i < data.size()) {
    		VariantDatum datum = data.get(i);
    		if(!dataIndex.containsKey(datum.loc.getContig())) {
    			Map<Integer, ArrayList<Integer>> posIndex = new HashMap<Integer, ArrayList<Integer>>();
    			dataIndex.put(datum.loc.getContig(), posIndex);
    		}
    		Map<Integer, ArrayList<Integer>> posIndex = dataIndex.get(datum.loc.getContig());
    		if(!posIndex.containsKey(datum.loc.getStart())) {
    			ArrayList<Integer> index = new ArrayList<Integer>();
    			posIndex.put(datum.loc.getStart(), index);
    		}
    		posIndex.get(datum.loc.getStart()).add(i);
    		i ++;
    	}
    }
    
    public List<Tranche> getTranches() {
    	return tranches;
    }
}
