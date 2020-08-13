package org.bgi.flexlab.gaea.tools.callsv;

public class Score {
	
	private final long LZERO = -10000000000L;
	private final long LSMALL = LZERO>>1;
	private final double minLogExp = -Math.log(-LZERO);
	private double logprob;
	
	public Score(){}
	
	public double logPoissionTailProb(int k, double lambda){
		logprob = LZERO;
		double plogprob;
		do{
			plogprob = logprob;
			logprob = LAdd(plogprob, logPoissionPDF(k++, lambda));
		}while(logprob - plogprob >= 0.01);
		
		return plogprob;
	}

	/*
	 * compute logP
	 */
	private double logPoissionPDF(int n, double lambda) {
		double logk_factorial = (n==0) ? 0 : (n*Math.log(n) - n + 0.5*Math.log(2*Math.PI*n));
		double log_lambda = (lambda <= 0) ? LSMALL : Math.log(lambda);
		double logp = n*log_lambda - lambda - logk_factorial;
		return logp ;
	}
	
	private double LAdd(double plogprob, double logp) {
		if(plogprob < logp){
			double tmp = plogprob;
			plogprob = logp;
			logp = tmp;
		}
		double diff = logp - plogprob;
		
		if(diff < minLogExp){
			return (plogprob < LSMALL) ? LZERO : plogprob;
		}else{
			return plogprob + Math.log(1+Math.pow(Math.E, diff));
		}
	}

}
