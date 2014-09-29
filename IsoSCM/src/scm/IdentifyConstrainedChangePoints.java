package scm;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import multisample.JointSegmentationResult;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;

import util.Util;
import util.Util.ExtremeObjectTracker;
import util.Util.IntegralTable;
import util.Util.MapList;
import cern.jet.random.NegativeBinomial;
import cern.jet.random.engine.RandomEngine;

public class IdentifyConstrainedChangePoints {



	/*
	 * Calculate the probability of every sequence from (s,t) of observations
	 * y coming from a distribution with a single parameter, integrating
	 * over possible parameter values with prior alpha and beta. 
	 */
	public static SCMResult calculate_segment_probabilities(double[] y, double alpha_0, double beta_0, double r){	

		IntegralTable cumulative = new IntegralTable(y);
		double[] log_binomial_y = new double[y.length];
		for(int i=0; i<y.length; i++){
			log_binomial_y[i] = Gamma.logGamma(y[i]+r+1)-(Gamma.logGamma(y[i]+1)+Gamma.logGamma(r+1));
		}
		IntegralTable logBinomial = new IntegralTable(log_binomial_y);

		double[][] log_p = new double[y.length][y.length];
		double[][] segment_mle = new double[y.length][y.length];

		double logBeta_alpha0_beta0 = Beta.logBeta(alpha_0, beta_0);

		for(int s=0; s<y.length; s++){
			for(int t=s; t<y.length; t++){
				int n = t-s+1;

				// compute the log marginal log likelihood
				double sum = cumulative.sum(s, t);
				double mean = sum/n;

				double alpha_1 = sum+alpha_0;
				double beta_1 = (n*r)+beta_0;

				double mle_p = (mean)/(r+mean);

				log_p[s][t] = Beta.logBeta(alpha_1, beta_1)-logBeta_alpha0_beta0+logBinomial.sum(s,t);
				segment_mle[s][t] = (mle_p*r)/(1-mle_p);
			}
		}

		return new SCMResult(log_p,segment_mle);
	}


	static class SCMResult {
		double[][] log_p;
		double[][] segment_mle;

		public SCMResult(double[][] log_p, double[][] segment_mle) {
			this.log_p = log_p;
			this.segment_mle = segment_mle;
		}
	}

	enum Regime {increasing, decreasing};

	public static SegmentationResult v_constrained_segmentation(double[] y, double alpha_0, double beta_0, double nb_r, int r, double p, double min_fold){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");

		double[] log_q	= new double[y.length];
		@SuppressWarnings("unchecked")
		// holds the likelihood for the segmentation starting at t, given cp at t-1
		ExtremeObjectTracker<Integer, Double>[] inc_map = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Integer, Double>[] dec_map = new ExtremeObjectTracker[y.length];
		// holds the extreme value
		ExtremeObjectTracker<Double, Double>[] inc_map_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] dec_map_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] inc_nxt_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] dec_nxt_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Regime, Double>[] nxt_regime = new ExtremeObjectTracker[y.length];
		
		SCMResult result = calculate_segment_probabilities(y, alpha_0, beta_0, nb_r);	

		double[][] log_p = result.log_p;
		double[][] segment_mle = result.segment_mle;

		int n = y.length;

		// length probability density function
		double[] log_g = new double[y.length];
		// length cumulative density function
		double[] log_G = new double[y.length];
		double[] g0 = new double[y.length];
		double[] G0 = new double[y.length];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<y.length; i++){
			inc_map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			inc_map_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			inc_nxt_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			dec_map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			dec_map_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			dec_nxt_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			nxt_regime[i] = new ExtremeObjectTracker<Regime, Double>(new Util.ComparableComparator<Double>());
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new LinkedList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];

				double cur = segment_mle[t][s];

				// if current > min of increasing
				// switching regimes from increasing to decreasing
				if(cur > inc_map_mle[s+1].getMaxObject()){
					double next_mle = inc_nxt_mle[s+1].getMaxObject();
					double fold = next_mle / cur;
					if(fold < min_fold){
					// if the fold is high enough... TODO
					double P = log_p[t][s]+inc_map[s+1].getMax()+ p_segment_length;
					dec_map[t].put(s+1, P);
					dec_map_mle[t].put(cur, P);
					dec_nxt_mle[t].put(segment_mle[t][s], P);
					nxt_regime[t].put(Regime.increasing, P);
					}
				}

				//maintain a decreasing regime
				if(dec_map_mle[s+1].hasExtrema() && cur > dec_map_mle[s+1].getMaxObject()){
					double next_mle = dec_nxt_mle[s+1].getMaxObject();
					double fold = next_mle / cur;
					if(fold < min_fold){
					// if the fold is high enough... TODO
					double P = log_p[t][s]+dec_map[s+1].getMax()+ p_segment_length;
					dec_map[t].put(s+1, P);
					dec_map_mle[t].put(cur, P);
					dec_nxt_mle[t].put(segment_mle[t][s], P);
					nxt_regime[t].put(Regime.decreasing, P);
					}
				}

				// maintain an increasing regime
				if(cur < inc_map_mle[s+1].getMaxObject()){
					
					double next_mle = inc_nxt_mle[s+1].getMaxObject();
					double fold = cur / next_mle;
					if(fold < min_fold){
					// if the fold is high enough... TODO
					double P = log_p[t][s]+inc_map[s+1].getMax()+ p_segment_length;
					inc_map[t].put(s+1, P);
					inc_map_mle[t].put(cur, P);
					inc_nxt_mle[t].put(segment_mle[t][s], P);
					nxt_regime[t].put(Regime.increasing, P);
					}
				}
//				double extreme = constrain_decreasing ? Math.max(inc_map_mle[s+1].getMaxObject(),segment_mle[t][s]) : Math.min(inc_map_mle[s+1].getMaxObject(),segment_mle[t][s]);

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t]);
			}

			// the probability of there being no more changepoints
			double P = log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1]));
			inc_map[t].put(n, P);
			inc_map_mle[t].put(segment_mle[t][n-1], P);
			inc_nxt_mle[t].put(segment_mle[t][n-1], P);
			nxt_regime[t].put(Regime.increasing, P);
			log_probs.add(P);

			//			map[t].put(n, log_p[t][n-1] + log_g[n-t-1]);
			//			map_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + log_g[n-t-1]);
			//			nxt_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + log_g[n-t-1]);
			//			log_probs.add(log_p[t][n-1] + log_g[n-t-1]);
			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();


		int first = 0;
		//		int last = inc_map[0].getExtreme();
		int last = -1;
		
	
		if(dec_map[0].hasExtrema() && dec_map[0].getMax() > inc_map[0].getMax())
			last = dec_map[0].getMaxObject();
		else
			last = inc_map[0].getMaxObject();
		
		Regime regime = nxt_regime[0].getExtreme();
		
		
		mapMLEs.add(segment_mle[first][last-1]);

		while(last<n){
//			System.out.printf("%d %d %s\n",first,last,regime);
			mapCPs.add(last-1);
			first = last;
			
			switch (regime) {
			case increasing:
				last = inc_map[last].getMaxObject();
				break;
			case decreasing:
				last = dec_map[last].getMaxObject();
				break;
			}
			
			regime = nxt_regime[first].getMaxObject();
			
			
			
			mapMLEs.add(segment_mle[first][last-1]);
		}


		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		return segmentation;
	}

	public static SegmentationResult a_constrained_segmentation(double[] y, double alpha_0, double beta_0, double nb_r, int r, double p, double min_fold){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");

		double[] log_q	= new double[y.length];
		@SuppressWarnings("unchecked")
		// holds the likelihood for the segmentation starting at t, given cp at t-1
		ExtremeObjectTracker<Integer, Double>[] inc_map = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Integer, Double>[] dec_map = new ExtremeObjectTracker[y.length];
		// holds the extreme value
		ExtremeObjectTracker<Double, Double>[] inc_map_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] dec_map_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] inc_nxt_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] dec_nxt_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Regime, Double>[] nxt_regime = new ExtremeObjectTracker[y.length];
		
		SCMResult result = calculate_segment_probabilities(y, alpha_0, beta_0, nb_r);	

		double[][] log_p = result.log_p;
		double[][] segment_mle = result.segment_mle;

		int n = y.length;

		// length probability density function
		double[] log_g = new double[y.length];
		// length cumulative density function
		double[] log_G = new double[y.length];
		double[] g0 = new double[y.length];
		double[] G0 = new double[y.length];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<y.length; i++){
			inc_map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			inc_map_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			inc_nxt_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			dec_map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			dec_map_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			dec_nxt_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			nxt_regime[i] = new ExtremeObjectTracker<Regime, Double>(new Util.ComparableComparator<Double>());
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new LinkedList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];

				double cur = segment_mle[t][s];

				// switching regimes from increasing to decreasing
				if(cur < dec_map_mle[s+1].getMaxObject()){
					double next_mle = dec_nxt_mle[s+1].getMaxObject();
					double fold = cur / next_mle;
					if(fold < min_fold){
					// if the fold is high enough... TODO
					double P = log_p[t][s]+dec_map[s+1].getMax()+ p_segment_length;
					inc_map[t].put(s+1, P);
					inc_map_mle[t].put(cur, P);
					inc_nxt_mle[t].put(segment_mle[t][s], P);
					nxt_regime[t].put(Regime.decreasing, P);
					}
				}

				//maintain a decreasing regime
				if(cur > dec_map_mle[s+1].getMaxObject()){
					double next_mle = dec_nxt_mle[s+1].getMaxObject();
					double fold = next_mle / cur;
					if(fold < min_fold){
					// if the fold is high enough... TODO
					double P = log_p[t][s]+dec_map[s+1].getMax()+ p_segment_length;
					dec_map[t].put(s+1, P);
					dec_map_mle[t].put(cur, P);
					dec_nxt_mle[t].put(segment_mle[t][s], P);
					nxt_regime[t].put(Regime.decreasing, P);
					}
				}

				// maintain an increasing regime
				if(inc_map_mle[s+1].hasExtrema() && cur < inc_map_mle[s+1].getMaxObject()){
					
					double next_mle = inc_nxt_mle[s+1].getMaxObject();
					double fold = cur / next_mle;
					if(fold < min_fold){
					// if the fold is high enough... TODO
					double P = log_p[t][s]+inc_map[s+1].getMax()+ p_segment_length;
					inc_map[t].put(s+1, P);
					inc_map_mle[t].put(cur, P);
					inc_nxt_mle[t].put(segment_mle[t][s], P);
					nxt_regime[t].put(Regime.increasing, P);
					}
				}

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t]);
			}

			// the probability of there being no more changepoints
			double P = log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1]));
			dec_map[t].put(n, P);
			dec_map_mle[t].put(segment_mle[t][n-1], P);
			dec_nxt_mle[t].put(segment_mle[t][n-1], P);
			nxt_regime[t].put(Regime.decreasing, P);
			log_probs.add(P);

			//			map[t].put(n, log_p[t][n-1] + log_g[n-t-1]);
			//			map_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + log_g[n-t-1]);
			//			nxt_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + log_g[n-t-1]);
			//			log_probs.add(log_p[t][n-1] + log_g[n-t-1]);
			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();


		int first = 0;
		//		int last = inc_map[0].getExtreme();
		int last = -1;
		
	
		if(inc_map[0].hasExtrema() && inc_map[0].getMax() > dec_map[0].getMax())
			last = inc_map[0].getMaxObject();
		else
			last = dec_map[0].getMaxObject();
		
		Regime regime = nxt_regime[0].getExtreme();
		
		
		mapMLEs.add(segment_mle[first][last-1]);

		while(last<n){
//			System.out.printf("%d %d %s\n",first,last,regime);
			mapCPs.add(last-1);
			first = last;
			
			switch (regime) {
			case increasing:
				last = inc_map[last].getMaxObject();
				break;
			case decreasing:
				last = dec_map[last].getMaxObject();
				break;
			}
			
			regime = nxt_regime[first].getMaxObject();
			
			
			
			mapMLEs.add(segment_mle[first][last-1]);
		}


		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		return segmentation;
	}
	
	/**
	 * 
	 * @param y
	 * @param alpha_0 - beta prior
	 * @param beta_0 - beta prior
	 * @param nb_r - neg binomial r parameter
	 * @param r - segment length r paramter
	 * @param p - segment length p parameter
	 * @return
	 */
	public static SegmentationResult constrained_viterbi_segmentation(double[] y, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");

		double[] log_q	= new double[y.length];
		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[] map = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] map_mle = new ExtremeObjectTracker[y.length];

		SCMResult result = calculate_segment_probabilities(y, alpha_0, beta_0, nb_r);

		double[][] log_p = result.log_p;
		double[][] segment_mle = result.segment_mle;

		int n = y.length;

		// length probability density function
		double[] log_g = new double[y.length];
		// length cumulative density function
		double[] log_G = new double[y.length];
		double[] g0 = new double[y.length];
		double[] G0 = new double[y.length];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<y.length; i++){
			map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			map_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new LinkedList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				//				if(t==0 || t==299)
				//				System.out.printf("(%d,%d) all %.1f < %.1f map s+1 %.1f %d\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject(),map[s+1].getMax(),map[s+1].getMaxObject());
				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])
				//					System.out.printf("(%d,%d) fail %.1f < %.1f\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject());

				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				double extreme = constrain_decreasing ? Math.max(map_mle[s+1].getMaxObject(),segment_mle[t][s]) : Math.min(map_mle[s+1].getMaxObject(),segment_mle[t][s]);
				//				if(map_mle[s+1].getMaxObject() > (segment_mle[t][s]+1) || map[s+1].getMax() == Double.NEGATIVE_INFINITY){

				if(constrain_decreasing){
					if(map_mle[s+1].getMaxObject() > (segment_mle[t][s])){
						map[t].put(s+1, Double.NEGATIVE_INFINITY);			
						map_mle[t].put(extreme, Double.NEGATIVE_INFINITY);
					}
					else{ 
						if(t>0){
							map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+log_g[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getExtremeValue()+log_g[s-t+1]);
						}
						else{
							map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+log_G[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getExtremeValue()+log_G[s-t+1]);	
						}
					}
				}
				// otherwise it should be constrained increasing
				else{
					if(map_mle[s+1].getMaxObject() < (segment_mle[t][s])){
						map[t].put(s+1, Double.NEGATIVE_INFINITY);			
						map_mle[t].put(extreme, Double.NEGATIVE_INFINITY);
					}
					else{ 
						if(t>0){
							map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+log_g[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getExtremeValue()+log_g[s-t+1]);
						}
						else{
							map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+log_G[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getExtremeValue()+log_G[s-t+1]);

						}
					}
				}

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t+1]);
			}

			// the probability of there being no more changepoints

			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			map_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_probs.add(log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();
		int first = 0;
		int last = map[0].getExtreme();
		mapMLEs.add(segment_mle[first][last-1]);
		//		System.out.println(map_mle[0].getMax());
		//		System.out.println(map_mle[0].getMaxObject());
		//		System.out.printf("last %d ll %.1f max %.1f\n", first, map_mle[first].getMax(), map_mle[first].getMaxObject());

		while(last<n){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));

			//			System.out.printf("last %d ll %.1f max %.1f\n", last, map_mle[last].getMax(), map_mle[last].getMaxObject());
			//			System.out.println(map_mle[last].getMax());
			//			System.out.println(map_mle[last].getMaxObject());
			mapCPs.add(last-1);
			first = last;
			last = map[last].getExtreme();
			mapMLEs.add(segment_mle[first][last-1]);
		}

		//		System.out.println(mapCPs);
		//		System.out.println(mapMLEs);

		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		//		System.out.println(Util.list(segmentation.segments()));
		//		mapCP.add(last);
		//		ensembl.add(mapCPs);

		//
		//		List<List<Integer>> ensembl = sample_ensembl(log_p, log_q, log_g, log_G, 50);
		//		plotEnsembl2(ensembl, 0, y.length-1);
		//		plotEnsembl3(ensembl, y);
		//			
		//		ChartUtils.showChart(ChartUtils.createLineChart(Util.map(Util.list("","1","2","4"),Util.enlist(data,data2,data3,data4)), "", "", "", false));
		return segmentation;
	}

	public static SegmentationResult doubly_constrained_segmentation(double[] y, double alpha_0, double beta_0, double nb_r, int r, double p, boolean constrain_decreasing, double min_fold){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");

		double[] log_q	= new double[y.length];
		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[] map = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] map_mle = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] nxt_mle = new ExtremeObjectTracker[y.length];

		SCMResult result = calculate_segment_probabilities(y, alpha_0, beta_0, nb_r);

		double[][] log_p = result.log_p;
		double[][] segment_mle = result.segment_mle;

		int n = y.length;

		// length probability density function
		double[] log_g = new double[y.length];
		// length cumulative density function
		double[] log_G = new double[y.length];
		double[] g0 = new double[y.length];
		double[] G0 = new double[y.length];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<y.length; i++){
			map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			map_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			nxt_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		//		double aBeta = 1;
		//		double bBeta = 8;
		//		double zBeta = Beta.logBeta(aBeta, bBeta);

		//TODO put log_G instead of log_g when t=0
		//TODO make sure log_G = logsum, not log_G += logsum
		//TODO make sure log_g = nb(i) not nb(i+1) so that sums to 1
		//TODO make sure that the length r is bigger than 1 to assign non-zero probs to long segments

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new LinkedList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				//				if(t==0 || t==299)
				//				System.out.printf("(%d,%d) all %.1f < %.1f map s+1 %.1f %d\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject(),map[s+1].getMax(),map[s+1].getMaxObject());
				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])
				//					System.out.printf("(%d,%d) fail %.1f < %.1f\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject());

				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				double next_mle = nxt_mle[s+1].getMaxObject();
				double extreme = constrain_decreasing ? Math.max(map_mle[s+1].getMaxObject(),segment_mle[t][s]) : Math.min(map_mle[s+1].getMaxObject(),segment_mle[t][s]);
				//				if(map_mle[s+1].getMaxObject() > (segment_mle[t][s]+1) || map[s+1].getMax() == Double.NEGATIVE_INFINITY){
				//				System.out.printf("(%d,%d) this %.2f next %.2f\n", t, s, segment_mle[t][s] , next_mle);

				Boolean fails_constraint = null;
				if(constrain_decreasing){
					double fold = next_mle/segment_mle[t][s];
					fails_constraint = map_mle[s+1].getMaxObject() >= (segment_mle[t][s]) || fold > min_fold;
				}
				else{
					double fold = segment_mle[t][s]/next_mle;
					fails_constraint = map_mle[s+1].getMaxObject() <= (segment_mle[t][s]) || fold > min_fold;
				}


				if(fails_constraint){
					map[t].put(s+1, Double.NEGATIVE_INFINITY);			
					map_mle[t].put(extreme, Double.NEGATIVE_INFINITY);
					nxt_mle[t].put(segment_mle[t][s], Double.NEGATIVE_INFINITY);
				}
				else{ 
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+ p_segment_length);
					map_mle[t].put(extreme, log_p[t][s]+map[s+1].getExtremeValue()+ p_segment_length);
					nxt_mle[t].put(segment_mle[t][s], log_p[t][s]+map[s+1].getExtremeValue()+ p_segment_length);
				}

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t]);
			}

			// the probability of there being no more changepoints

			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			map_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			nxt_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_probs.add(log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			
//			map[t].put(n, log_p[t][n-1] + log_g[n-t-1]);
//			map_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + log_g[n-t-1]);
//			nxt_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + log_g[n-t-1]);
//			log_probs.add(log_p[t][n-1] + log_g[n-t-1]);
			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();
		int first = 0;
		int last = map[0].getExtreme();
		mapMLEs.add(segment_mle[first][last-1]);
		//		System.out.println(map_mle[0].getMax());
		//		System.out.println(map_mle[0].getMaxObject());
		//		System.out.printf("last %d ll %.1f max %.1f\n", first, map_mle[first].getMax(), map_mle[first].getMaxObject());

		while(last<n){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));

			//			System.out.printf("last %d ll %.1f max %.1f\n", last, map_mle[last].getMax(), map_mle[last].getMaxObject());
			//			System.out.println(map_mle[last].getMax());
			//			System.out.println(map_mle[last].getMaxObject());
			mapCPs.add(last-1);
			first = last;
			last = map[last].getExtreme();
			mapMLEs.add(segment_mle[first][last-1]);
		}

		/*
		 * calculate the increase in likelihood for adding a change point
		if(mapCPs.size()>0){
			
		for(int i=0; i<mapCPs.size(); i++){
			int prev = -1;
			int next = -1;
			if(i==0){
				prev=0;
			}
			else{
				prev = mapCPs.get(i-1);
			}
			if(i==mapCPs.size()-1){
				next = n-1;
			}
			else{
				next = mapCPs.get(i+1);
			}
			System.out.printf("prev %d next %d i %d p %.3e\n", prev,next,i,(log_p[prev+1][mapCPs.get(i)]+log_p[mapCPs.get(i)+1][next]-log_p[prev+1][next]));
			
			System.out.println(i+" "+(log_p[prev+1][mapCPs.get(i)]+log_p[mapCPs.get(i)+1][next]-log_p[prev+1][next]));
		}
		}
		*/
		//		System.out.println(mapCPs);
		//		System.out.println(mapMLEs);

		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		//		System.out.println(Util.list(segmentation.segments()));
		//		mapCP.add(last);
		//		ensembl.add(mapCPs);

		//
		//		List<List<Integer>> ensembl = sample_ensembl(log_p, log_q, log_g, log_G, 50);
		//		plotEnsembl2(ensembl, 0, y.length-1);
		//		plotEnsembl3(ensembl, y);
		//			
		//		ChartUtils.showChart(ChartUtils.createLineChart(Util.map(Util.list("","1","2","4"),Util.enlist(data,data2,data3,data4)), "", "", "", false));
		return segmentation;
	}

	public static JointSegmentationResult doubly_constrained_multi_segmentation(double[][] y, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing, double min_fold){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");
		int n = y[0].length;

		double[] log_q	= new double[n];
		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[] map = new ExtremeObjectTracker[n];
		ExtremeObjectTracker<Double, Double>[][] map_mle = new ExtremeObjectTracker[y.length][n];
		ExtremeObjectTracker<Double, Double>[][] nxt_mle = new ExtremeObjectTracker[y.length][n];
		double[][][] log_p = new double[y.length][][];
		double[][][] segment_mle = new double[y.length][][];

		SCMResult[] results = new SCMResult[y.length];
		for(int i=0; i<y.length; i++){
			results[i] = calculate_segment_probabilities(y[i], alpha_0, beta_0, nb_r);
			log_p[i] = results[i].log_p;
			segment_mle[i] = results[i].segment_mle;
		}




		// length probability density function
		double[] log_g = new double[n];
		// length cumulative density function
		double[] log_G = new double[n];
		double[] g0 = new double[n];
		double[] G0 = new double[n];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<n; i++){
			map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			for(int j=0; j<y.length; j++){
				map_mle[j][i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
				nxt_mle[j][i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			}
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		//		double aBeta = 1;
		//		double bBeta = 8;
		//		double zBeta = Beta.logBeta(aBeta, bBeta);

		//TODO put log_G instead of log_g when t=0
		//TODO make sure log_G = logsum, not log_G += logsum
		//TODO make sure log_g = nb(i) not nb(i+1) so that sums to 1
		//TODO make sure that the length r is bigger than 1 to assign non-zero probs to long segments

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new LinkedList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				//				if(t==0 || t==299)
				//				System.out.printf("(%d,%d) all %.1f < %.1f map s+1 %.1f %d\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject(),map[s+1].getMax(),map[s+1].getMaxObject());
				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])
				//					System.out.printf("(%d,%d) fail %.1f < %.1f\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject());

				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				boolean all_decreasing = true;
				boolean any_satisfy_fold = false;

				double log_p_t_s = 0;

				for(int i=0; i<y.length; i++){
					log_p_t_s += log_p[i][t][s];

					double next_mle = nxt_mle[i][s+1].getMaxObject();
					double extreme = constrain_decreasing ? Math.max(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]) : Math.min(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]);
					//				if(map_mle[s+1].getMaxObject() > (segment_mle[t][s]+1) || map[s+1].getMax() == Double.NEGATIVE_INFINITY){
					//				System.out.printf("(%d,%d) this %.2f next %.2f\n", t, s, segment_mle[t][s] , next_mle);

					
					Boolean satisfies_decreasing = null;
					Boolean satisfies_fold = null;
					if(constrain_decreasing){
						double fold = next_mle/segment_mle[i][t][s];
						satisfies_fold = fold <= min_fold;
						satisfies_decreasing = map_mle[i][s+1].getMaxObject() < (segment_mle[i][t][s]);
					}
					else{
						double fold = segment_mle[i][t][s]/next_mle;
						satisfies_fold = fold <= min_fold;
						satisfies_decreasing = map_mle[i][s+1].getMaxObject() > (segment_mle[i][t][s]);
					}
					all_decreasing &= satisfies_decreasing;
					any_satisfy_fold |= satisfies_fold;
				}

				for(int i=0; i<y.length; i++){
					double extreme = constrain_decreasing ? Math.max(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]) : Math.min(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]);
						if(!(all_decreasing&&any_satisfy_fold)){
						map_mle[i][t].put(extreme, Double.NEGATIVE_INFINITY);
						nxt_mle[i][t].put(segment_mle[i][t][s], Double.NEGATIVE_INFINITY);
					}
					else{ 
						double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
						map_mle[i][t].put(extreme, log_p_t_s+map[s+1].getExtremeValue()+ p_segment_length);
						nxt_mle[i][t].put(segment_mle[i][t][s], log_p_t_s+map[s+1].getExtremeValue()+ p_segment_length);
					}
				}

				if(!(all_decreasing&&any_satisfy_fold)){
					map[t].put(s+1, Double.NEGATIVE_INFINITY);
				}
				else{
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[t].put(s+1, log_p_t_s+map[s+1].getExtremeValue()+ p_segment_length);
				}


				// sum over the positions of the next change point
				log_probs.add(log_p_t_s+log_q[s+1]+log_g[s-t]);
			}

			double log_p_t_s = 0;

			for(int i=0; i<y.length; i++){
				log_p_t_s += log_p[i][t][n-1];
			}
			// the probability of there being no more changepoints

			map[t].put(n, log_p_t_s + Math.log(1 - Math.exp(log_G[n-t-1])));
			for(int i=0; i<y.length; i++){
			map_mle[i][t].put(segment_mle[i][t][n-1], log_p_t_s + Math.log(1 - Math.exp(log_G[n-t-1])));
			nxt_mle[i][t].put(segment_mle[i][t][n-1], log_p_t_s + Math.log(1 - Math.exp(log_G[n-t-1])));
			}
			log_probs.add(log_p_t_s + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double[]> mapMLEs = new LinkedList<Double[]>();
		List<Double> mapLs = new LinkedList<Double>();
		int first = 0;
		int last = map[0].getExtreme();
		
		Double[] mle = new Double[y.length];
		for(int i=0; i<y.length; i++){
			mle[i] = segment_mle[i][first][last-1];
		}
		
		mapMLEs.add(mle);
		//		System.out.println(map_mle[0].getMax());
		//		System.out.println(map_mle[0].getMaxObject());
		//		System.out.printf("last %d ll %.1f max %.1f\n", first, map_mle[first].getMax(), map_mle[first].getMaxObject());

		while(last<n){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));

			//			System.out.printf("last %d ll %.1f max %.1f\n", last, map_mle[last].getMax(), map_mle[last].getMaxObject());
			//			System.out.println(map_mle[last].getMax());
			//			System.out.println(map_mle[last].getMaxObject());
			mapCPs.add(last-1);
			first = last;
			last = map[last].getExtreme();
			mle = new Double[y.length];
			
			for(int i=0; i<y.length; i++){
				mle[i] = segment_mle[i][first][last-1];
			}
			
			mapMLEs.add(mle);
		}

		//		System.out.println(mapCPs);
		//		System.out.println(mapMLEs);
		
		double[][] mles = new double[mapMLEs.size()][y.length];
		for(int i=0; i<mapMLEs.size(); i++){
		for(int j=0; j<y.length; j++){
			mles[i][j] = mapMLEs.get(i)[j];
		}
		}
		

		JointSegmentationResult segmentation = new JointSegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])),	mles, ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		//		System.out.println(Util.list(segmentation.segments()));
		//		mapCP.add(last);
		//		ensembl.add(mapCPs);

		//
		//		List<List<Integer>> ensembl = sample_ensembl(log_p, log_q, log_g, log_G, 50);
		//		plotEnsembl2(ensembl, 0, y.length-1);
		//		plotEnsembl3(ensembl, y);
		//			
		//		ChartUtils.showChart(ChartUtils.createLineChart(Util.map(Util.list("","1","2","4"),Util.enlist(data,data2,data3,data4)), "", "", "", false));
		return segmentation;
	}

	public static SegmentationResult fold_constrained_segmentation(double[] y, double alpha_0, double beta_0, int nb_r, int r, double p, double min_fold){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");

		double[] log_q	= new double[y.length];
		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[] map = new ExtremeObjectTracker[y.length];
		ExtremeObjectTracker<Double, Double>[] nxt_mle = new ExtremeObjectTracker[y.length];

		SCMResult result = calculate_segment_probabilities(y, alpha_0, beta_0, nb_r);

		double[][] log_p = result.log_p;
		double[][] segment_mle = result.segment_mle;

		int n = y.length;

		// length probability density function
		double[] log_g = new double[y.length];
		// length cumulative density function
		double[] log_G = new double[y.length];
		double[] g0 = new double[y.length];
		double[] G0 = new double[y.length];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<y.length; i++){
			map[i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			nxt_mle[i] = new ExtremeObjectTracker<Double, Double>(new Util.ComparableComparator<Double>());
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		//		double aBeta = 1;
		//		double bBeta = 8;
		//		double zBeta = Beta.logBeta(aBeta, bBeta);
		//		double min_fold = .5;

		//TODO put log_G instead of log_g when t=0
		//TODO make sure log_G = logsum, not log_G += logsum
		//TODO make sure log_g = nb(i) not nb(i+1) so that sums to 1
		//TODO make sure that the length r is bigger than 1 to assign non-zero probs to long segments

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new LinkedList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				//				if(t==0 || t==299)
				//				System.out.printf("(%d,%d) all %.1f < %.1f map s+1 %.1f %d\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject(),map[s+1].getMax(),map[s+1].getMaxObject());
				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])
				//					System.out.printf("(%d,%d) fail %.1f < %.1f\n", t,s, segment_mle[t][s], map_mle[s+1].getMaxObject());

				//				if(map_mle[s+1].getMaxObject() > segment_mle[t][s])

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				double next_mle = nxt_mle[s+1].getMaxObject();

				Boolean fails_constraint = null;

				if(Math.max(next_mle,segment_mle[t][s])==0){
					fails_constraint=true;
				}
				else{
					double fold = Math.min(next_mle,segment_mle[t][s])/Math.max(next_mle,segment_mle[t][s]);
					fails_constraint = fold > min_fold;
				}



				if(fails_constraint){
					map[t].put(s+1, Double.NEGATIVE_INFINITY);			
					nxt_mle[t].put(segment_mle[t][s], Double.NEGATIVE_INFINITY);
				}
				else{ 
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+ p_segment_length);
					nxt_mle[t].put(segment_mle[t][s], log_p[t][s]+map[s+1].getExtremeValue()+ p_segment_length);
				}

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t]);
			}

			// the probability of there being no more changepoints

			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			nxt_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_probs.add(log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();
		int first = 0;
		int last = map[0].getExtreme();
		mapMLEs.add(segment_mle[first][last-1]);

		while(last<n){
			mapCPs.add(last-1);
			first = last;
			last = map[last].getExtreme();
			mapMLEs.add(segment_mle[first][last-1]);
		}

		//		System.out.println(mapCPs);
		//		System.out.println(mapMLEs);

		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		//		System.out.println(Util.list(segmentation.segments()));
		//		mapCP.add(last);
		//		ensembl.add(mapCPs);

		//
		//		List<List<Integer>> ensembl = sample_ensembl(log_p, log_q, log_g, log_G, 50);
		//		plotEnsembl2(ensembl, 0, y.length-1);
		//		plotEnsembl3(ensembl, y);
		//			
		//		ChartUtils.showChart(ChartUtils.createLineChart(Util.map(Util.list("","1","2","4"),Util.enlist(data,data2,data3,data4)), "", "", "", false));
		return segmentation;
	}

	/**
	 * 
	 * @param y
	 * @param alpha_0
	 * @param beta_0
	 * @param nb_r
	 * @param r
	 * @param p
	 * @param m the number of change points
	 * @return
	 */
	public static SegmentationResult viterbi_segmentation(double[] y, double alpha_0, double beta_0, int nb_r, int r, double p, int m){
		//		m=m+1;
		double[][] log_q_j	= new double[m+1][y.length];
		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[][] map_j = new ExtremeObjectTracker[m+1][y.length];

		SCMResult result = calculate_segment_probabilities(y, alpha_0, beta_0, nb_r);
		double[][] log_p = result.log_p;
		double[][] segment_mle = result.segment_mle;

		int n = y.length;

		// length probability density function
		double[] log_g = new double[y.length];
		// length cumulative density function
		double[] log_G = new double[y.length];
		double[] g0 = new double[y.length];
		double[] G0 = new double[y.length];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<y.length; i++){
			for(int j=0; j<m+1; j++){
				map_j[j][i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			}
			log_g[i] = log_nb(i+1, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] += Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		// perform the recursions
		{
			int j=m+1-1;
			for(int t = n-1; t > j-1; t--){
				int s = n-1;
				log_q_j[j][t] = log_p[t][s]+Math.log(1 - Math.exp(log_G[s-t+1]));
				map_j[j][t].put(s+1, log_p[t][s]+Math.log(1 - Math.exp(log_G[s-t+1])));
			}

		}
		for(int j=m+1-2; j>-1; j--){
			for(int t = n-m-1+j; t > j-1; t--){
				List<Double> log_probs = new LinkedList<Double>();
				for(int s=t; s<n-m-1+j+1; s++){
					log_probs.add(log_p[t][s]+log_q_j[j+1][s+1]+log_g[s-t+1]);
					map_j[j][t].put(s+1, log_p[t][s]+map_j[j+1][s+1].getExtremeValue()+log_g[s-t+1]);
				}
				log_q_j[j][t] = Util.logSum(log_probs);
			}
		}

		/*
		for(int t=n-1; t>-1; t--){

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point
				map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+log_g[s-t+1]);
				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q_j[s+1]+log_g[s-t+1]);
			}

			// the probability of there being no more changepoints
			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_probs.add(log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_q_j[t] = Util.logSum(log_probs);
		}
		 */

		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();
		int first = 0;
		int j=0;
		int last = map_j[j][0].getExtreme();
		mapMLEs.add(segment_mle[first][last-1]);
		j++;

		while(j<m+1){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));
			mapCPs.add(last-1);
			first = last;
			last = map_j[j][last].getExtreme();
			mapMLEs.add(segment_mle[first][last-1]);
			j++;
		}

		//		System.out.println(mapCPs);
		//		System.out.println(mapMLEs);

		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		//		System.out.println(Util.list(segmentation.segments()));
		//		mapCP.add(last);
		//		ensembl.add(mapCPs);

		//
		//		List<List<Integer>> ensembl = sample_ensembl(log_p, log_q, log_g, log_G, 50);
		//		plotEnsembl2(ensembl, 0, y.length-1);
		//		plotEnsembl3(ensembl, y);
		//			
		//		ChartUtils.showChart(ChartUtils.createLineChart(Util.map(Util.list("","1","2","4"),Util.enlist(data,data2,data3,data4)), "", "", "", false));
		return segmentation;
	}

	public static SegmentationResult viterbi_segmentation(double[][] y, double alpha_0, double beta_0, int nb_r, int r, double p, int m){
		//		m=m+1;

		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[][] map_j = new ExtremeObjectTracker[m+1][y[0].length];

		SCMResult[] results = new SCMResult[y.length];
		for(int i=0; i<y.length; i++){
			results[i] = calculate_segment_probabilities(y[i], alpha_0, beta_0, nb_r);
		}

		//		(y, alpha_0, beta_0, nb_r);
		//		double[][] log_p = result.log_p;
		//		double[][] segment_mle = result.segment_mle;

		int n = y[0].length;

		double[][] log_q_j	= new double[m+1][n];
		// length probability density function
		double[] log_g = new double[n];
		// length cumulative density function
		double[] log_G = new double[n];
		double[] g0 = new double[n];
		double[] G0 = new double[n];

		//TODO can integrate over the length of the first segment like Fearnhead

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<n; i++){
			for(int j=0; j<m+1; j++){
				map_j[j][i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			}
			log_g[i] = log_nb(i+1, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] += Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		// perform the recursions
		{
			int j=m+1-1;
			for(int t = n-1; t > j-1; t--){
				int s = n-1;
				double log_p_t_s = 0;
				for(int i=0; i<results.length; i++){
					log_p_t_s += results[i].log_p[t][s];

				}
				log_q_j[j][t] = log_p_t_s+Math.log(1 - Math.exp(log_G[s-t+1]));
				map_j[j][t].put(s+1, log_p_t_s+Math.log(1 - Math.exp(log_G[s-t+1])));
			}	
		}
		for(int j=m+1-2; j>-1; j--){
			for(int t = n-m-1+j; t > j-1; t--){
				List<Double> log_probs = new LinkedList<Double>();
				for(int s=t; s<n-m-1+j+1; s++){
					double log_p_t_s = 0;
					for(int i=0; i<results.length; i++){
						log_p_t_s += results[i].log_p[t][s];
					}
					log_probs.add(log_p_t_s+log_q_j[j+1][s+1]+log_g[s-t+1]);
					map_j[j][t].put(s+1, log_p_t_s+map_j[j+1][s+1].getExtremeValue()+log_g[s-t+1]);
				}
				log_q_j[j][t] = Util.logSum(log_probs);
			}
		}

		/*
		for(int t=n-1; t>-1; t--){

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point
				map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+log_g[s-t+1]);
				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q_j[s+1]+log_g[s-t+1]);
			}

			// the probability of there being no more changepoints
			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_probs.add(log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_q_j[t] = Util.logSum(log_probs);
		}
		 */

		//		double d = Math.exp(log_p[previousChangePoint+1][i]+log_q[i+1]+log_g[i-previousChangePoint] 
		//				- log_q[previousChangePoint+1]);			
		//		states.add(i);
		//		probabilities.add(d);

		MapList<String, Double[]> dat = new MapList<String, Double[]>();

		try{
			for(int i=0; i<n; i++){
				int previousChangePoint = -1;
				//				double d = Math.exp(results[0].log_p[previousChangePoint+1][i]+log_q_j[0][i+1]+log_g[i-previousChangePoint] 
				//						- log_q_j[0][previousChangePoint+1]);	
				double d = (results[0].log_p[previousChangePoint+1][i]+results[1].log_p[previousChangePoint+1][i]+results[2].log_p[previousChangePoint+1][i]+log_q_j[1][i+1]+log_g[i-previousChangePoint] 
						- log_q_j[0][previousChangePoint+1]);	

				dat.put("1", new Double[]{i+0.,d});
				System.out.println(d);
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		try{
			int previousChangePoint = 68;
			for(int i=previousChangePoint; i<n; i++){
				//				double d = Math.exp(results[0].log_p[previousChangePoint+1][i]+log_q_j[0][i+1]+log_g[i-previousChangePoint] 
				//						- log_q_j[0][previousChangePoint+1]);	
				double d = (results[0].log_p[previousChangePoint+1][i]+results[1].log_p[previousChangePoint+1][i]+results[2].log_p[previousChangePoint+1][i]+log_q_j[2][i+1]+log_g[i-previousChangePoint] 
						- log_q_j[1][previousChangePoint+1]);	

				dat.put("2", new Double[]{i+0.,d});
				System.out.println(d);
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		try{
			int previousChangePoint = 241;
			for(int i=previousChangePoint; i<n; i++){
				//				double d = Math.exp(results[0].log_p[previousChangePoint+1][i]+log_q_j[0][i+1]+log_g[i-previousChangePoint] 
				//						- log_q_j[0][previousChangePoint+1]);	
				double d = (results[0].log_p[previousChangePoint+1][i]+results[1].log_p[previousChangePoint+1][i]+results[2].log_p[previousChangePoint+1][i]+log_q_j[3][i+1]+log_g[i-previousChangePoint] 
						- log_q_j[2][previousChangePoint+1]);	

				dat.put("3", new Double[]{i+0.,d});
				System.out.println(d);
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		
		//		List<List<Integer>> ensembl = new LinkedList<List<Integer>>();
		List<Integer> mapCPs = new LinkedList<Integer>();
		List<Double> mapMLEs = new LinkedList<Double>();
		List<Double> mapLs = new LinkedList<Double>();
		int first = 0;
		int j=0;
		int last = map_j[j][0].getExtreme();
		//		mapMLEs.add(segment_mle[first][last-1]);
		j++;

		while(j<m+1){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));
			mapCPs.add(last-1);
			first = last;
			last = map_j[j][last].getExtreme();
			//			mapMLEs.add(segment_mle[first][last-1]);
			j++;
		}

		//		System.out.println(mapCPs);
		//		System.out.println(mapMLEs);

		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);
		//		System.out.println(Util.list(segmentation.segments()));
		//		mapCP.add(last);
		//		ensembl.add(mapCPs);

		//
		//		List<List<Integer>> ensembl = sample_ensembl(log_p, log_q, log_g, log_G, 50);
		//		plotEnsembl2(ensembl, 0, y.length-1);
		//		plotEnsembl3(ensembl, y);
		//			
		//		ChartUtils.showChart(ChartUtils.createLineChart(Util.map(Util.list("","1","2","4"),Util.enlist(data,data2,data3,data4)), "", "", "", false));
		return segmentation;
	}


	/**
	 * 
	 * Probability of getting x successes and before r failures
	 * when the probability of success is p in each Bernoulli
	 * trial
	 * 
	 * @param x - number of successes
	 * @param r - number of failures
	 * @param p - probability of success
	 * @return
	 */
	public static double log_nb(int x, int r, double p){
		//		MathUtils.binomialCoefficient(r+x, x);

		return Gamma.logGamma(r+x) - Gamma.logGamma(x+1) - Gamma.logGamma(r) + ((r)*Math.log(1-p)) + (x*Math.log(p));
	}
	
	public static double log_nb_real(int x, double r, double p){
		//		MathUtils.binomialCoefficient(r+x, x);

		return Gamma.logGamma(r+x) - Gamma.logGamma(x+1) - Gamma.logGamma(r) + ((r)*Math.log(1-p)) + (x*Math.log(p));
	}

	public static void main(String[] args) throws IOException {
//		testRealSmoothing();
//		testRealCtn5();
//		testMultiConstrained();
		//		System.out.println(Math.log(1-Math.exp(1e-5)));
		//		System.out.printf("%e\n", 1-Math.exp(1e-5));
		//		System.exit(0);
		//	testRealShank2();
		//		testRealCars();
//		testRealSerinc3();
		//		testFlyData();
		//				testRealData();

		int l=50;
		double[][] Y2 = new double[2][l*3];
		double[] Y = new double[l*4];
		MapList<String, Double[]> data = new MapList<String, Double[]>();
		//		RandomEngine rng =  new DRand();
		RandomEngine rng = RandomEngine.makeDefault();

		List<Double> Ys = new LinkedList<Double>();


		int R=5;
		NegativeBinomial pd = new NegativeBinomial(R, .51, rng);

		//		NormalDistributionImpl nd = new NormalDistributionImpl(0, 1);
		for(int i=0; i<l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, .75, rng);
		for(int i=l; i<2*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, .77, rng);
		for(int i=2*l; i<3*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, .65, rng);
		for(int i=3*l; i<4*l; i++){
			Y[i] = pd.nextDouble();
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		l=(3*l)/2;
		pd = new NegativeBinomial(10, .15, rng);
		for(int i=0; i<1*l; i++){
			Y2[1][i] = pd.nextDouble();
			//			data.put("2", new Double[]{(double) i,Y2[1][i]});
		}
		pd = new NegativeBinomial(100, .01, rng);
		for(int i=1*l; i<2*l; i++){
			Y2[1][i] = pd.nextDouble();
			//			data.put("2", new Double[]{(double) i,Y2[1][i]});
		}

		pd = new NegativeBinomial(R, .75, rng);
		for(int i=0; i<50; i++){
			Ys.add(pd.nextDouble());
		}
		pd = new NegativeBinomial(R, .45, rng);
		for(int i=0; i<10; i++){
			Ys.add(pd.nextDouble());
		}
		pd = new NegativeBinomial(R, .8, rng);
		for(int i=0; i<50; i++){
			Ys.add(pd.nextDouble());
		}
		pd = new NegativeBinomial(R, .1, rng);
		for(int i=0; i<10; i++){
			Ys.add(pd.nextDouble());
		}

		data.getMap().clear();
		{
			double i=0;
			for(Double d : Ys){
				data.put("1", new Double[]{i,d});
				i++;
			}
			Y = ArrayUtils.toPrimitive(Ys.toArray(new Double[0]));
		}
		//TODO set mu_0 as mean, and expected variance as overall variance?
		double m_0=5;
		double V_0=1;
		double a_0=1;
		double b_0=.001;
		//		System.out.println(Util.list(Y));
		//		for(int i=0; i<3*l-1; i++){
		//			System.out.println(Y.length);
		//			System.out.println(i+1);
		//			System.out.println(2*l);
		//			double[] y1=Arrays.copyOfRange(Y, 0, i+1);
		//			double[] y2=Arrays.copyOfRange(Y, i+1, 3*l);
		//			double ll1 = log_marginal_likelihood(m_0, V_0, a_0, b_0, y1);
		//			double ll2 = log_marginal_likelihood(m_0, V_0, a_0, b_0, y2);
		//			System.out.println(Util.list(y1)+":"+Util.list(y2));
		//			System.out.println(ll1);
		//			System.out.println(ll2);
		//			data.put("ll", new Double[]{(double) i,ll1+ll2});
		//		}


		//		DynamicSegmentProbabilities dsp = new DynamicSegmentProbabilities(Y, m_0, V_0, a_0, b_0, 1, .95);
		double alpha_0=1; 
		double beta_0=1; 
		int nb_r = 30;

		boolean constrained_decreasing = true;
		//		SegmentationResult sr = fold_constrained_segmentation(Y, alpha_0, beta_0, nb_r, 2, .99, .9);
				SegmentationResult sr = doubly_constrained_segmentation(Y, alpha_0, beta_0, nb_r, 50, .99, constrained_decreasing, .9);
		//		SegmentationResult sr = viterbi_segmentation(Y, alpha_0, beta_0, nb_r, 1, .95, 3);
//		SegmentationResult sr = BetaNegativeBinomial.viterbi_segmentation(Y, alpha_0, beta_0, R, 25, .95);
		//		SegmentationResult sr = viterbi_segmentation(Y, m_0, V_0, a_0, b_0, 1, .95);
		//		System.out.println(Util.list(sr.change_point));
		MapList<String,Double[]> cps = new MapList<String, Double[]>();
		for(int i=0; i<sr.change_point.length; i++){
			//			data.put("a", new Double[]{sr.change_point[i]+0.,0.});
			data.put("a", new Double[]{sr.change_point[i]+0.,-1.});
			cps.put("", new Double[]{sr.change_point[i]+0.,0.});
			cps.put("", new Double[]{sr.change_point[i]+0.,1.});
			cps.put("", new Double[]{sr.change_point[i]+0.,null});
		}

//		JFreeChart jfc0 = ChartUtils.createScatterChart(data.getMap(), "", "", "", false);
//		ChartUtils.setDomainAxisRange(jfc0, 0, Ys.size());
//		ChartUtils.showChart(jfc0);
//		//		ChartUtils.savePDF(jfc0, "/home/sol/workspace/Latex/isoscm/figures/unconstrained_data.pdf", 300, 200);
//		JFreeChart jfc1 = ChartUtils.createLineChart(cps.getMap(), "", "", "", false);
//		ChartUtils.setDomainAxisRange(jfc1, 0, Ys.size());
//		ChartUtils.showChart(jfc1);
		//		ChartUtils.savePDF(jfc1, "/home/sol/workspace/Latex/isoscm/figures/unconstrained_cp.pdf", 300, 100);

		//		System.out.println(dynamic_marginal_loglikelihood(Util.list(sr.change_point), dsp));
		//	
		//		for(int i=0; i<dsp.l-1; i++){
		//			System.out.printf("%d %f\n", i, dynamic_marginal_loglikelihood(Util.list(i), dsp));
		//		}
		//		System.out.printf("%f\n", dynamic_marginal_loglikelihood(Util.list(), dsp));
		//		
		////		System.exit(0);
		//		
		//		for(int i=0; i<dsp.l-1; i++){
		//			for(int j=i+1; j<dsp.l-1; j++){
		//					System.out.printf("%d %d %f\n", i,j, dynamic_marginal_loglikelihood(Util.list(i,j), dsp));
		//			}
		//		}

	}
}
