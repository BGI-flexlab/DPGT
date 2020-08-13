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
import htsjdk.variant.variantcontext.VariantContextUtils;

import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocationParser;
import org.bgi.flexlab.gaea.tools.mapreduce.vcfqualitycontrol.VCFQualityControlOptions;
import org.bgi.flexlab.gaea.util.MathUtils;
import org.bgi.flexlab.gaea.util.QualityUtils;
import org.bgi.flexlab.gaea.util.RandomUtils;

import java.io.IOException;

public class VariantDatumMessenger{

	private double[] annotations;
	private boolean[] isNull;
	private double lod;
	
	/**
	 * isKnown:atTruthSite:atTrainingSite:atAntiTrainingSite:isTransition:isSNP:failingSTDThreshold
	 */
	public static final byte isKnown = 0x1;
	public static final byte atTruthSite = 0x1 << 1;
	public static final byte atTrainingSite = 0x1 << 2;
	public static final byte atAntiTrainingSite = 0x1 << 3;
	public static final byte isTransition = 0x1 << 4;
	public static final byte isSNP = 0x1 << 5;
	public static final byte failingSTDThreshold = 0x1 << 6;
	private byte flag;

	private double originalQual;
	private double prior;
	private int consensusCount;
	private GenomeLocation loc;
	private int worstAnnotation;
	
	private VariantDatumMessenger(Builder builder) {
		this.annotations = builder.annotations;
		this.isNull = builder.isNull;
		this.lod = builder.lod;
		this.flag = builder.flag;
		this.originalQual = builder.originalQual;
		this.prior = builder.prior;
		this.consensusCount = builder.consensusCount;
		this.loc = builder.loc;
		this.worstAnnotation = builder.worstAnnotation;
	}
	
	public String toString() {
		StringBuilder datumInfo = new StringBuilder();
		datumInfo.append(annotations.length);
		for(int i = 0; i < isNull.length; i++) {
			datumInfo.append("\t");
			datumInfo.append(isNull[i]);
		}
		
		for(int i = 0; i < annotations.length; i++) {
			datumInfo.append("\t");
			datumInfo.append(annotations[i]);
		}
		
		datumInfo.append("\t");datumInfo.append(lod);
		datumInfo.append("\t");datumInfo.append(flag);
		datumInfo.append("\t");datumInfo.append(originalQual);
		datumInfo.append("\t");datumInfo.append(prior);
		datumInfo.append("\t");datumInfo.append(consensusCount);
		datumInfo.append("\t");datumInfo.append(loc.getContig());
		datumInfo.append("\t");datumInfo.append(loc.getStart());
		datumInfo.append("\t");datumInfo.append(loc.getStop());
		datumInfo.append("\t");datumInfo.append(worstAnnotation);
		
		return datumInfo.toString();
	}
	
	public double[] getAnnotations() {
		return annotations;
	}
	
	public boolean[] getIsNull() {
		return isNull;
	}
	
	public double getLod() {
		return lod;
	}
	
	public byte getFlag() {
		return flag;
	}
	
	public boolean checkFlag(byte i) {
		return (flag & i) != 0;
	}
	  
	public static boolean checkFlag(byte flag, byte i) {
		return (flag & i) != 0;
	}
	public double getOriginalQual() {
		return originalQual;
	}
	
	public double getPrior() {
		return prior;
	}
	
	public int getConsensusCount() {
		return consensusCount;
	}
	
	public GenomeLocation getLoc() {
		return loc;
	}
	
	public int getWorstAnnotation() {
		return worstAnnotation;
	}
	
	public static class Builder {
		
		private double[] annotations;
		private boolean[] isNull;
		private double lod;
	
		private byte flag;

		private double originalQual;
		private double prior;
		private int consensusCount;
		private GenomeLocation loc;
		private int worstAnnotation;
		
		private VCFQualityControlOptions options;
		
		private ResourceManager manager;
				
		private VariantContext vc;
				
//		builder used for retrieving object
		public Builder(){
			this.annotations = null;
			this.isNull = null;
			this.lod = 0.0;
			this.flag = 0;
			this.originalQual = 0;
			this.prior = 2.0;
			this.consensusCount = 0;
			this.loc = null;
			this.worstAnnotation = 0;
		}
		
		public Builder(ResourceManager manager, VariantContext vc, VCFQualityControlOptions options) {
		// TODO Auto-generated constructor stub
			this();
			this.vc = vc;
			this.options = options;
			this.manager = manager;
		}
		
		public Builder decodeAnnotations() {
			decodeAnnotations(vc, true);
			return this;
		}
		
		public Builder setFlagV() throws IOException {
			if(vc.isSNP() && vc.isBiallelic())
				this.flag = (this.flag |= VariantDatumMessenger.isSNP);
			if(checkFlag(this.flag, VariantDatumMessenger.isSNP) && VariantContextUtils.isTransition(vc))
				this.flag = (this.flag |= VariantDatumMessenger.isTransition);

			this.flag = parseTrainingSets(flag, options.isTrustAllPolymorphic());
			return this;
		}
		
		public Builder setOriginalQual() {
			this.originalQual = vc.getPhredScaledQual();
			return this;
		}
		
		public Builder setPrior() {
			double priorFactor = QualityUtils.qualityToProbability(this.prior);
			this.prior = Math.log10( priorFactor ) - Math.log10( 1.0 - priorFactor );
			return this;
		}
		
		public Builder setLoc(GenomeLocationParser genomeLocParser) {
			this.loc = genomeLocParser.createGenomeLocation(vc);
			return this;
		}
		
		private void decodeAnnotations( final VariantContext vc, final boolean jitter) {
	        int iii = 0;
	        isNull = new boolean[options.getUseAnnotations().size()];
	        annotations = new double[options.getUseAnnotations().size()];
	        for( final String key : options.getUseAnnotations() ) {
	            isNull[iii] = false;
	            annotations[iii] = decodeAnnotation( key, vc, jitter );
	            if( Double.isNaN(annotations[iii]) ) { isNull[iii] = true; }
	            iii++;
	        }
	    }

	    private static double decodeAnnotation( final String annotationKey, final VariantContext vc, final boolean jitter ) {
	        double value;

	        try {
	            value = vc.getAttributeAsDouble( annotationKey, Double.NaN );
	            if( Double.isInfinite(value) ) { value = Double.NaN; }
	            if( jitter && annotationKey.equalsIgnoreCase("HRUN") ) { // Integer valued annotations must be jittered a bit to work in this GMM
	                  value += -0.25 + 0.5 * RandomUtils.getRandomGenerator().nextDouble();
	            }

	            if( jitter && annotationKey.equalsIgnoreCase("HaplotypeScore") && MathUtils.compareDoubles(value, 0.0, 0.0001) == 0 ) { value = -0.2 + 0.4*RandomUtils.getRandomGenerator().nextDouble(); }
	            if( jitter && annotationKey.equalsIgnoreCase("FS") && MathUtils.compareDoubles(value, 0.0, 0.001) == 0 ) { value = -0.2 + 0.4*RandomUtils.getRandomGenerator().nextDouble(); }
	        } catch( Exception e ) {
	            value = Double.NaN; // The VQSR works with missing data by marginalizing over the missing dimension when evaluating the Gaussian mixture model
	        }

	        return value;
	    }
	    
	    private byte parseTrainingSets(byte flag, final boolean TRUST_ALL_POLYMORPHIC ) throws IOException {
	        byte result = flag;
	    	for( final TrainData trainingData : manager.getTrainDataSet() ) {
	            for( final VariantContext trainVC : trainingData.get(loc) ) {
	            	if( isValidVariant( vc, trainVC, TRUST_ALL_POLYMORPHIC ) ) {
	                    if(checkFlag(this.flag, VariantDatumMessenger.isKnown) || trainingData.isKnown()) {
	                    	result = (result |= VariantDatumMessenger.isKnown);
	                    }
	                    if(checkFlag(this.flag, VariantDatumMessenger.atTruthSite) || trainingData.isTruth()) {
	                    	result = (result |= VariantDatumMessenger.atTruthSite);
	                    }
	                    if(checkFlag(this.flag, VariantDatumMessenger.atTrainingSite) || trainingData.isTraining()) {
	                    	result = (result |= VariantDatumMessenger.atTrainingSite);
	                    }
	                    this.prior = Math.max(this.prior, trainingData.getPrior());
	                }
	                if( trainVC != null ) {
	                	if(checkFlag(this.flag, VariantDatumMessenger.atAntiTrainingSite) || trainingData.isAntiTraining()) {
	                    	result = (result |= VariantDatumMessenger.atAntiTrainingSite);
	                    }
	                }
	            }
	        }
	    	return result;
	    }
	    
	    private boolean isValidVariant( final VariantContext evalVC, final VariantContext trainVC, final boolean TRUST_ALL_POLYMORPHIC) {
	        return trainVC != null && trainVC.isNotFiltered() && trainVC.isVariant() && ResourceManager.checkVariationClass( evalVC, trainVC ) &&
	                        (TRUST_ALL_POLYMORPHIC || !trainVC.hasGenotypes() || trainVC.isPolymorphicInSamples());
	    }
	    
	    public VariantDatumMessenger buildFrom(String line, GenomeLocationParser genomeLocParser) {
			String[] lineSplits = line.split("\t");
			
			int index = 0;
			//int id = Integer.parseInt(lineSplits[index++]);
			int annoSize = Integer.parseInt(lineSplits[index++]);
			this.isNull = new boolean[annoSize];
			for(int i = 0; i < annoSize; i++) {
				this.isNull[i] = Boolean.parseBoolean(lineSplits[index++]);
			}
			
			this.annotations = new double[annoSize];
			for(int i = 0; i < annoSize; i++) {
				if(!this.isNull[i])
					this.annotations[i] = Double.parseDouble(lineSplits[index++]);
				else
					index++;
			}
			
			this.lod = Double.parseDouble(lineSplits[index++]);
			this.flag = Byte.parseByte(lineSplits[index++]);
			this.originalQual = Double.parseDouble(lineSplits[index++]);
			this.prior = Double.parseDouble(lineSplits[index++]);
			this.consensusCount = Integer.parseInt(lineSplits[index++]);
			this.loc = genomeLocParser.createGenomeLocation(lineSplits[index++], Integer.parseInt(lineSplits[index++]), Integer.parseInt(lineSplits[index++]));
			this.worstAnnotation = Integer.parseInt(lineSplits[index++]);
			//return id;
			return new VariantDatumMessenger(this);
		}
	    
		public VariantDatumMessenger build() {
			return new VariantDatumMessenger(this);
		}
	}
	
}
