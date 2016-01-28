package scm;

import java.util.ArrayList;
import java.util.List;

import util.Util;

public class ConfidenceCalculator {

	double[] log_q_prefix;
	double[] log_q_suffix;
	
	List<Double> likelihoods;
	
	public ConfidenceCalculator(PrefixSuffixResult psr) {
		this.log_q_prefix = psr.log_q_prefix;
		this.log_q_suffix = psr.log_q_suffix;
		likelihoods = new ArrayList<Double>();
		for(int i=0; i<psr.l-1; i++){
			likelihoods.add(log_q_prefix[i]+log_q_suffix[i+1]);
		}
	}
	
	/**
	 * 
	 * @param start
	 * @param end
	 * @return confidence that observation at position start and end are on different segments (the confidence that there is a change point at some position between them)
	 */
	public ConfidenceResult calculateConfidence(int start, int end){
		double log_probility = Util.logSum(likelihoods.subList(start, end))-log_q_suffix[0];
		double log_odds = Math.exp(0-log_probility)-1 <=0 ? Double.POSITIVE_INFINITY : log_probility - Util.logSubtract(0, log_probility); 
		return new ConfidenceResult(log_probility, log_odds);
	}
	
	public ConfidenceResult calculateConfidence(int p){
		double log_probility = likelihoods.get(p)-log_q_suffix[0];
		double log_odds = Math.exp(0-log_probility)-1 <=0 ? Double.POSITIVE_INFINITY : log_probility - Util.logSubtract(0, log_probility); 
		return new ConfidenceResult(log_probility, log_odds);
	}
}
