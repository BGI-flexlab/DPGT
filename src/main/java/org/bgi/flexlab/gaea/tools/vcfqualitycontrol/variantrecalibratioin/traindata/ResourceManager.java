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
package org.bgi.flexlab.gaea.tools.vcfqualitycontrol.variantrecalibratioin.traindata;

import htsjdk.variant.variantcontext.VariantContext;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;

import java.util.ArrayList;
import java.util.List;

public class ResourceManager {
	private final double[] annotations;
    private final boolean[] isNull;
    public final List<String> annotationKeys;
    //protected final List<TrainingSet> trainingSets;
    protected final List<TrainData> trainingSets;
    
    public ResourceManager(VCFQualityControlOptions options) {
        this.annotationKeys = new ArrayList<String>( options.getUseAnnotations() );
        annotations = new double[annotationKeys.size()];
        isNull = new boolean[annotationKeys.size()];
        trainingSets = new ArrayList<TrainData>();
    }

    public void addTrainingSet( final TrainData trainingSet ) {
        trainingSets.add( trainingSet );
    }

    public boolean checkHasTrainingSet() {
        for( final TrainData trainingSet : trainingSets ) {
            if( trainingSet.isTraining() ) { return true; }
        }
        return false;
    }

    public boolean checkHasTruthSet() {
        for( final TrainData trainingSet : trainingSets ) {
            if( trainingSet.isTruth() ) { return true; }
        }
        return false;
    }

    public boolean checkHasKnownSet() {
        for( final TrainData trainingSet : trainingSets ) {
            if( trainingSet.isKnown() ) { return true; }
        }
        return false;
    }
    
    public double[] getAnnotations() {
    	return annotations;
    }
    
    public boolean[] getIsNull() {
    	return isNull;
    }
    
    public List<TrainData> getTrainDataSet() {
    	return trainingSets;
    }

    public static boolean checkVariationClass( final VariantContext evalVC, final VariantContext trainVC ) {
        switch( trainVC.getType() ) {
            case SNP:
            case MNP:
                return checkVariationClass( evalVC, VCFQualityControlOptions.Mode.SNP );
            case INDEL:
            case MIXED:
            case SYMBOLIC:
                return checkVariationClass( evalVC, VCFQualityControlOptions.Mode.INDEL );
            default:
                return false;
        }
    }

    public static boolean checkVariationClass( final VariantContext evalVC, final VCFQualityControlOptions.Mode mode ) {
        switch( mode ) {
            case SNP:
                return evalVC.isSNP() || evalVC.isMNP();
            case INDEL:
                return evalVC.isStructuralIndel() || evalVC.isIndel() || evalVC.isMixed() || evalVC.isSymbolic();
            case BOTH:
                return true;
            default:
                throw new RuntimeException( "Encountered unknown recal mode: " + mode );
        }
    }

}
