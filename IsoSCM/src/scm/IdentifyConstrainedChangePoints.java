package scm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import multisample.JointSegmentationResult;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jfree.chart.JFreeChart;

import util.ChartUtils;
import util.Util;
import util.Util.ExtremeObjectTracker;
import util.Util.IntegralTable;
import util.Util.MapCounter;
import util.Util.MapList;
import cern.jet.random.NegativeBinomial;
import cern.jet.random.engine.RandomEngine;

public class IdentifyConstrainedChangePoints {

	final static Logger l = LogManager.getLogger();

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
			List<Double> log_probs = new ArrayList<Double>();

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

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();


		int first = 0;
		//		int last = inc_map[0].getMaxObject();
		int last = -1;


		if(dec_map[0].hasExtrema() && dec_map[0].getMax() > inc_map[0].getMax())
			last = dec_map[0].getMaxObject();
		else
			last = inc_map[0].getMaxObject();

		Regime regime = nxt_regime[0].getMaxObject();


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
			List<Double> log_probs = new ArrayList<Double>();

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

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();


		int first = 0;
		//		int last = inc_map[0].getExtreme();
		int last = -1;


		if(inc_map[0].hasExtrema() && inc_map[0].getMax() > dec_map[0].getMax())
			last = inc_map[0].getMaxObject();
		else
			last = dec_map[0].getMaxObject();

		Regime regime = nxt_regime[0].getMaxObject();


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
			List<Double> log_probs = new ArrayList<Double>();

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
							map[t].put(s+1, log_p[t][s]+map[s+1].getMaxObject()+log_g[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getMaxObject()+log_g[s-t+1]);
						}
						else{
							map[t].put(s+1, log_p[t][s]+map[s+1].getMaxObject()+log_G[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getMaxObject()+log_G[s-t+1]);	
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
							map[t].put(s+1, log_p[t][s]+map[s+1].getMaxObject()+log_g[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getMaxObject()+log_g[s-t+1]);
						}
						else{
							map[t].put(s+1, log_p[t][s]+map[s+1].getMaxObject()+log_G[s-t+1]);
							map_mle[t].put(extreme, log_p[t][s]+map[s+1].getMaxObject()+log_G[s-t+1]);

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

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();
		int first = 0;
		int last = map[0].getMaxObject();
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
			last = map[last].getMaxObject();
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

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new ArrayList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				double next_mle = nxt_mle[s+1].getMaxObject();
				double extreme = constrain_decreasing ? Math.max(map_mle[s+1].getMaxObject(),segment_mle[t][s]) : Math.min(map_mle[s+1].getMaxObject(),segment_mle[t][s]);

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
					map[t].put(s+1, log_p[t][s]+map[s+1].getMax()+ p_segment_length);
					map_mle[t].put(extreme, log_p[t][s]+map[s+1].getMax()+ p_segment_length);
					nxt_mle[t].put(segment_mle[t][s], log_p[t][s]+map[s+1].getMax()+ p_segment_length);
				}

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t]);
			}

			// the probability of there being no more changepoints

			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			map_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			nxt_mle[t].put(segment_mle[t][n-1], log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
			log_probs.add(log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));

			log_q[t] = Util.logSum(log_probs);
		}

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();
		int first = 0;
		int last = map[0].getMaxObject();
		mapMLEs.add(segment_mle[first][last-1]);

		while(last<n){
			mapCPs.add(last-1);
			first = last;
			last = map[last].getMaxObject();
			mapMLEs.add(segment_mle[first][last-1]);
		}


		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);

		return segmentation;
	}

	public static PrefixSuffixResult prefix_suffix(double[] y, double alpha_0, double beta_0, int nb_r, int r, double p){
		return prefix_suffix(new double[][]{y}, alpha_0, beta_0, nb_r, r, p);
	}

	public static PrefixSuffixResult prefix_suffix(double[][] y, double alpha_0, double beta_0, int nb_r, int r, double p){
		int l = y[0].length;

		double[][] log_p = new double[l][l];
		double[][][] segment_mle = new double[y.length][][];

		SCMResult[] results = new SCMResult[y.length];
		for(int i=0; i<y.length; i++){
			results[i] = calculate_segment_probabilities(y[i], alpha_0, beta_0, nb_r);
			for(int j=0; j<l; j++){
				for(int k=0; k<l; k++){
					log_p[j][k] += results[i].log_p[j][k];
				}
			}

			segment_mle[i] = results[i].segment_mle;
		}

		// length probability density function
		double[] log_g = new double[l];
		// length cumulative density function
		double[] log_G = new double[l];

		log_G[0] = Double.NEGATIVE_INFINITY;

		// pre-compute the length functions
		for(int i=0; i<l; i++){
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
		}

		double[] log_q_suffix	= new double[l];



		// perform the recursions
		for(int t=l-1; t>-1; t--){
			List<Double> log_probs = new ArrayList<Double>();

			for(int s=t; s<l-1; s++){
				// find the most likely position of the positions of the next change point

				double log_p_t_s = log_p[t][s];

				double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];	

				log_probs.add(log_p_t_s+log_q_suffix[s+1]+p_segment_length);
			}

			double log_p_t_s = log_p[t][l-1];

			double p_segment_length = Math.log(1 - Math.exp(log_G[l-t-1]));

			// the probability of there being no more changepoints
			log_probs.add(log_p_t_s + p_segment_length);

			log_q_suffix[t] = Util.logSum(log_probs);		
		}

		double[] log_q_prefix	= new double[l];

		for(int s=0; s<l; s++){

			List<Double> log_probs = new ArrayList<Double>();

			for(int t=s; t>0; t--){
				double log_p_t_s = log_p[t][s];

				double p_segment_length = s==l-1? Math.log(1 - Math.exp(log_G[l-t-1])) : log_g[s-t];	

				log_probs.add(log_p_t_s+log_q_prefix[t-1]+p_segment_length);
			}

			double log_p_t_s = log_p[0][s];

			double p_segment_length = Math.log(1 - Math.exp(log_G[s]));

			// the probability of there being no more changepoints
			log_probs.add(log_p_t_s + p_segment_length);

			log_q_prefix[s] = Util.logSum(log_probs);		
		}

		double[] log_suffix = new double[l];
		double[] log_prefix = new double[l];

		for(int t=l-1; t>-1; t--){
			double log_p_t_s = log_p[t][l-1];
			double p_segment_length = Math.log(1 - Math.exp(log_G[l-t-1]));
			log_suffix[t] = log_p_t_s + p_segment_length;		
		}

		for(int s=0; s<l; s++){
			double log_p_t_s = log_p[0][s];
			double p_segment_length = Math.log(1 - Math.exp(log_G[s]));
			log_prefix[s] = log_p_t_s + p_segment_length;		
		}

		return new PrefixSuffixResult(log_p, log_g, log_G, log_q_prefix, log_q_suffix, log_prefix, log_suffix, l);

	}

	public static int[] sample(PrefixSuffixResult psr, int n){
		int[] sample = suffix_sample(psr.log_p, psr.log_q_suffix, psr.log_g, psr.log_G, psr.l, n);
//		MapList<String, Integer[]> ml = new MapList<String, Integer[]>();
//		for (int i = 0; i < sample.length; i++) {
//			ml.put("", new Integer[]{i,sample[i]});
//		}
//		JFreeChart jfc = ChartUtils.createScatterChart(ml.getMap(), "", "", "", false);
//		ChartUtils.showChart(jfc);
		return sample;
	}

	public static int[] suffix_sample(double[][] log_p, double[] log_q, double[] log_g, double[] log_G, int l, int n){
		int[] n_with_cp = new int[l+1];
		n_with_cp[0] = n;

		for(int i=0; i<l; i++){

			int nt = n_with_cp[i];
			if(nt > 0){
				List<ChangePointSample> numNextChangePoints = sample_next_segment_start(log_p, log_q, log_g, log_G, i, nt, l);

				for(ChangePointSample cps : numNextChangePoints){
					n_with_cp[cps.segmentStart] += cps.n;
				}
			}
		}

		return n_with_cp;
	}

	static class ChangePointSample{
		// number of times segment starting at this index was sampled
		int n;
		// index of the first observation in the segment
		int segmentStart;
		public ChangePointSample(int n, int segmentStart) {
			this.n = n;
			this.segmentStart = segmentStart;
		}
	}
	/**
	 * 
	 * @param log_p
	 * @param log_q
	 * @param log_g
	 * @param log_G
	 * @param segmentStart
	 * @param n
	 * @param l
	 * @return an list of pairs of segment starts and how many times that change-point was sampled
	 */
	public static List<ChangePointSample> sample_next_segment_start(double[][] log_p, double[] log_q, double[] log_g, double[] log_G, int segmentStart, int n, int l){

		{
			int t = segmentStart;
			int[] next_segment_start = new int[l-t];
			double[] next_segment_start_probability = new double[l-t];

			for(int s=t; s<l-1; s++){
				// find the most likely position of the positions of the next change point

				double log_p_t_s = log_p[t][s];

				double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];	

				next_segment_start[s-t] = s+1;
				next_segment_start_probability[s-t] = Math.exp(log_p_t_s+log_q[s+1]+p_segment_length-log_q[t]);
			}

			double log_p_t_s = log_p[t][l-1];

			double p_segment_length = Math.log(1 - Math.exp(log_G[l-t-1]));

			// the probability of there being no more changepoints
			next_segment_start[l-1-t] = l;
			next_segment_start_probability[l-1-t] = Math.exp(log_p_t_s + p_segment_length-log_q[t]);

			//			System.out.printf("sum2 %f\n", Util.sum(Util.list(next_segment_start_probability)));
			int[] sample_n = sample_n(next_segment_start_probability, n);

			List<ChangePointSample> cps = new ArrayList<ChangePointSample>();
			for(int i=0; i<sample_n.length; i++){
				int nextSegmentStart = i+t+1;
				cps.add(new ChangePointSample(sample_n[i], nextSegmentStart));
			}

			return cps;

			/*
			// perform the recursions
			for(int t=l-1; t>-1; t--){
				List<Double> log_probs = new LinkedList<Double>();

				for(int s=t; s<l-1; s++){
					// find the most likely position of the positions of the next change point

					double log_p_t_s = log_p[t][s];

					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];	

					log_probs.add(log_p_t_s+log_q[s+1]+p_segment_length);
				}

				double log_p_t_s = log_p[t][l-1];

				double p_segment_length = Math.log(1 - Math.exp(log_G[l-t-1]));

				// the probability of there being no more changepoints
				log_probs.add(log_p_t_s + p_segment_length);

				log_q[t] = Util.logSum(log_probs);		
			}
			 */
		}

		//		double[] probabilities = new double[l];
		//
		//		// all possible next change points
		//		for(int i=segmentStart+1; i<l-1; i++){
		//
		//			double log_p_s_t =  log_p[segmentStart+1][i];
		//
		//			probabilities[i] = Math.exp(log_p_s_t+log_q[i+1]+log_g[i-segmentStart]-log_q[segmentStart+1]);
		//		}
		//
		//		// probability of no more change points
		//		{
		//
		//			double log_p_s_n = log_p[segmentStart+1][l-1];
		//			
		//			probabilities[l-1] = Math.exp(log_p_s_n-log_q[segmentStart+1])*(1-Math.exp(log_G[l-(segmentStart+1)-1]));
		//		}
		//		double[] primitive = ArrayUtils.toPrimitive(Util.normalize(Util.list(probabilities)).toArray(new Double[0]));
		//		//TODO probabilities must sum to 1
		//		System.out.println(Util.list(primitive));
		//		return sample_n(primitive, n);
	}

	/**
	 * 
	 * @param P_tau - pmf for a change point at time t
	 * @param n - number of change points to sample
	 */
	public static int[] sample_n(double[] P_tau, int n){

		int[] samples = new int[P_tau.length];

		ExponentialDistribution ed = new ExponentialDistribution(1);

		double[] x_i = new double[n+1];
		double[] s_i = new double[n+1];
		double S = 0;
		for (int i = 0; i < n+1; i++) {
			double x = ed.sample();
			x_i[i] = x;
			S+=x;
			s_i[i]=S;
		}
		double[] u_i = new double[n+1];
		for (int i = 0; i < n+1; i++) {
			u_i[i] = s_i[i]/S;
		}

		double Q=0;
		double U=u_i[0];
		int j=0,i=0;
		while(i<n){
			if(U<Q+P_tau[j]){
				samples[j]++;
				U = u_i[i+1];
				i++;
			}
			else{
				Q+=P_tau[j]; 
				j++;
			}
		}

		return samples;
	}

	public static JointSegmentationResult doubly_constrained_multi_segmentation(double[][] y, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing, double min_fold){
		int n = y[0].length;

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

		// perform the recursions
		for(int t=n-1; t>-1; t--){

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

				/*
				 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
				 * the next segment is 0, 
				 */

				boolean any_satisfy_constraints = false;
				boolean all_positive_upstream = true;

				double log_p_t_s = 0;

				double[] extremes = new double[y.length];
				for(int i=0; i<y.length; i++){
					log_p_t_s += log_p[i][t][s];

					double next_mle = nxt_mle[i][s+1].getMaxObject();
					extremes[i] = constrain_decreasing ? Math.max(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]) : Math.min(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]);
					//					
					Boolean satisfies_decreasing = null;
					Boolean satisfies_fold = null;

					if(constrain_decreasing){
						double fold = next_mle/segment_mle[i][t][s];
						satisfies_fold = fold <= min_fold;
						satisfies_decreasing = map_mle[i][s+1].getMaxObject() < (segment_mle[i][t][s]);
						all_positive_upstream &= segment_mle[i][t][s] > 0;
					}
					else{
						double fold = segment_mle[i][t][s]/next_mle;
						satisfies_fold = fold <= min_fold;
						satisfies_decreasing = (segment_mle[i][t][s]) < map_mle[i][s+1].getMaxObject();
						all_positive_upstream &= next_mle > 0;
					}
					any_satisfy_constraints |= (satisfies_fold && satisfies_decreasing);
				}

				if(!(any_satisfy_constraints&&all_positive_upstream)){
					map[t].put(s+1, Double.NEGATIVE_INFINITY);			
					for(int i=0; i<y.length; i++){
						map_mle[i][t].put(extremes[i], Double.NEGATIVE_INFINITY);
						nxt_mle[i][t].put(segment_mle[i][t][s], Double.NEGATIVE_INFINITY);
					}
				}
				else{
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[t].put(s+1, log_p_t_s+map[s+1].getMax()+ p_segment_length);
					for(int i=0; i<y.length; i++){
						map_mle[i][t].put(extremes[i], log_p_t_s+map[s+1].getMax()+ p_segment_length);
						nxt_mle[i][t].put(segment_mle[i][t][s], log_p_t_s+map[s+1].getMax()+ p_segment_length);
					}
				}

				/*
				 * 	if(fails_constraint){
					map[t].put(s+1, Double.NEGATIVE_INFINITY);			
					map_mle[t].put(extreme, Double.NEGATIVE_INFINITY);
					nxt_mle[t].put(segment_mle[t][s], Double.NEGATIVE_INFINITY);
				}
				else{ 
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[t].put(s+1, log_p[t][s]+map[s+1].getMaxObject()+ p_segment_length);
					map_mle[t].put(extreme, log_p[t][s]+map[s+1].getMaxObject()+ p_segment_length);
					nxt_mle[t].put(segment_mle[t][s], log_p[t][s]+map[s+1].getMaxObject()+ p_segment_length);
				}
				 */

				for(int i=0; i<y.length; i++){
					double extreme = constrain_decreasing ? Math.max(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]) : Math.min(map_mle[i][s+1].getMaxObject(),segment_mle[i][t][s]);
					if(!(any_satisfy_constraints&&all_positive_upstream)){
						map_mle[i][t].put(extreme, Double.NEGATIVE_INFINITY);
						nxt_mle[i][t].put(segment_mle[i][t][s], Double.NEGATIVE_INFINITY);
					}
					else{ 
						double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
						map_mle[i][t].put(extreme, log_p_t_s+map[s+1].getMax()+ p_segment_length);
						nxt_mle[i][t].put(segment_mle[i][t][s], log_p_t_s+map[s+1].getMax()+ p_segment_length);
					}
				}

				if(!(any_satisfy_constraints&&all_positive_upstream)){
					map[t].put(s+1, Double.NEGATIVE_INFINITY);
				}
				else{
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[t].put(s+1, log_p_t_s+map[s+1].getMax()+ p_segment_length);
				}

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
		}

		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double[]> mapMLEs = new ArrayList<Double[]>();
		List<Double> mapLs = new ArrayList<Double>();
		int first = 0;
		int last = map[0].getMaxObject();

		Double[] mle = new Double[y.length];
		for(int i=0; i<y.length; i++){
			mle[i] = segment_mle[i][first][last-1];
		}

		mapMLEs.add(mle);

		while(last<n){
			mapCPs.add(last-1);
			first = last;
			last = map[last].getMaxObject();
			mle = new Double[y.length];

			for(int i=0; i<y.length; i++){
				mle[i] = segment_mle[i][first][last-1];
			}

			mapMLEs.add(mle);
		}

		double[][] mles = new double[mapMLEs.size()][y.length];
		for(int i=0; i<mapMLEs.size(); i++){
			for(int j=0; j<y.length; j++){
				mles[i][j] = mapMLEs.get(i)[j];
			}
		}

		JointSegmentationResult segmentation = new JointSegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])),	mles, ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);

		return segmentation;
	}

	/**
	 * 
	 * @param nSegments
	 * @param y
	 * @param alpha_0
	 * @param beta_0
	 * @param nb_r
	 * @param r
	 * @param p
	 * @param constrain_decreasing
	 * @param min_fold
	 * @return change-point for a fixed number of segments
	 */
	public static JointSegmentationResult fixed_n_joint_segmentation(int nSegments, double[][] y, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing, double min_fold){
		//		if(!constrain_decreasing)
		//			throw new RuntimeException("unimplemented");
		int n = y[0].length;

		@SuppressWarnings("unchecked")
		// The most likely start position of the segment following a changepoint at position i
		ExtremeObjectTracker<Integer, Double>[][] map = new ExtremeObjectTracker[nSegments][n];
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
			log_g[i] = log_nb(i, r, p);
			if(i>0)
				log_G[i] = log_G[i-1];
			l.info("log_G[i]={}",log_G[i]);
			l.info("log_g[i]={}",log_g[i]);
			log_G[i] = Util.logSum(Util.list(log_G[i],log_g[i]));
			l.info("2 log_G[i]={}",log_g[i]);
		}
		// initialize datastructures
		for(int m=0; m<nSegments; m++){
			for(int i=0; i<n; i++){
				map[m][i] = new ExtremeObjectTracker<Integer, Double>(new Util.ComparableComparator<Double>());
			}
		}

		// calculate for last segment
		{
			int m=nSegments-1;
			for(int t=n-1; t>m-1; t--){
				l.info("{}th segment start t={}",m+2,t);

				double log_p_t_s = 0;
				double p_segment_length = Math.log(1 - Math.exp(log_G[n-t]));
				l.info("segment probability {}",p_segment_length);

				for(int i=0; i<y.length; i++){
					l.info("length log_p[i][t] = {}",log_p[i][t].length);
					l.info("likelihood sample {} = {}",i,log_p[i][t][n-1]);
					log_p_t_s += log_p[i][t][n-1];
				}

				double log_likelihood = p_segment_length + log_p_t_s;

				l.info("likelihood last segment starts at {} = {}",t+1,log_likelihood);

				map[m][t].put(n, log_likelihood);
			}
		}

		l.info("sequence length n={}",n);
		for(int m=nSegments-2; m>-1; m--){
			l.info("change point index m={}",m);
			// perform the recursions
			for(int t=n-1; t>m-1; t--){
				l.info("change point index m={}",m);
				l.info("segment start, the index of the previous change point t={}",t);

				l.info("s<{}",n-1-(nSegments-m));
				for(int s=t; s<n-1-(nSegments-m-1); s++){
					l.info("change point index m={}",m);
					l.info("segment start t={}",t);
					l.info("segment end s={}",s);

					// find the most likely position of the positions of the next change point

					/*
					 * The constraint: if the MLE of the next segment is greater than this segment, or the probability of 
					 * the next segment is 0, 
					 */

					boolean any_satisfy_fold = false;

					double log_p_t_s = 0;

					for(int i=0; i<y.length; i++){
						log_p_t_s += log_p[i][t][s];

					}

					//				for(int i=0; i<y.length; i++){
					//						if(!(any_satisfy_fold)){
					//						map_mle[m][i][t].put(extreme, Double.NEGATIVE_INFINITY);
					//						nxt_mle[m][i][t].put(segment_mle[i][t][s], Double.NEGATIVE_INFINITY);
					//					}
					//					else{ 
					//						double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					//					}
					//				}
					l.info("is null? => {}",map[m+1][s+1].getMax());

					//				if(!(any_satisfy_fold)){
					//					map[m][t].put(s+1, Double.NEGATIVE_INFINITY);
					//				}
					//				else{
					double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
					map[m][t].put(s+1, log_p_t_s+map[m+1][s+1].getMax()+ p_segment_length);
					//				}

				}

			}
		}

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double[]> mapMLEs = new ArrayList<Double[]>();
		List<Double> mapLs = new ArrayList<Double>();
		int m=0;
		int first = 0;
		int last = map[m][0].getMaxObject();
		l.info("traceback m={}",m);
		l.info("traceback max index for m last={}",map[m][0].getMaxObject());
		l.info("traceback max value={}",map[m][0].getMax());

		Double[] mle = new Double[y.length];
		for(int i=0; i<y.length; i++){
			mle[i] = segment_mle[i][first][last-1];
		}

		mapMLEs.add(mle);

		for(m=1; m<nSegments; m++){
			l.info("traceback m={}",m);
			l.info("traceback max index for m last={}",map[m][last].getMaxObject());
			l.info("traceback max value={}",map[m][last].getMax());
			mapCPs.add(last-1);
			first = last;
			last = map[m][last].getMaxObject();
			mle = new Double[y.length];

			for(int i=0; i<y.length; i++){
				mle[i] = segment_mle[i][first][last-1];
			}

			mapMLEs.add(mle);
		}

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

		// perform the recursions
		for(int t=n-1; t>-1; t--){
			List<Double> log_probs = new ArrayList<Double>();

			for(int s=t; s<n-1; s++){
				// find the most likely position of the positions of the next change point

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
					map[t].put(s+1, log_p[t][s]+map[s+1].getMax()+ p_segment_length);
					nxt_mle[t].put(segment_mle[t][s], log_p[t][s]+map[s+1].getMax()+ p_segment_length);
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

		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();
		int first = 0;
		int last = map[0].getMaxObject();
		mapMLEs.add(segment_mle[first][last-1]);

		while(last<n){
			mapCPs.add(last-1);
			first = last;
			last = map[last].getMaxObject();
			mapMLEs.add(segment_mle[first][last-1]);
		}

		SegmentationResult segmentation = new SegmentationResult(ArrayUtils.toPrimitive(mapCPs.toArray(new Integer[0])), ArrayUtils.toPrimitive(mapMLEs.toArray(new Double[0])), ArrayUtils.toPrimitive(mapLs.toArray(new Double[0])), n);

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
				List<Double> log_probs = new ArrayList<Double>();
				for(int s=t; s<n-m-1+j+1; s++){
					log_probs.add(log_p[t][s]+log_q_j[j+1][s+1]+log_g[s-t+1]);
					map_j[j][t].put(s+1, log_p[t][s]+map_j[j+1][s+1].getMaxObject()+log_g[s-t+1]);
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

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();
		int first = 0;
		int j=0;
		int last = map_j[j][0].getMaxObject();
		mapMLEs.add(segment_mle[first][last-1]);
		j++;

		while(j<m+1){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));
			mapCPs.add(last-1);
			first = last;
			last = map_j[j][last].getMaxObject();
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
				List<Double> log_probs = new ArrayList<Double>();
				for(int s=t; s<n-m-1+j+1; s++){
					double log_p_t_s = 0;
					for(int i=0; i<results.length; i++){
						log_p_t_s += results[i].log_p[t][s];
					}
					log_probs.add(log_p_t_s+log_q_j[j+1][s+1]+log_g[s-t+1]);
					map_j[j][t].put(s+1, log_p_t_s+map_j[j+1][s+1].getMaxObject()+log_g[s-t+1]);
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

		//		List<List<Integer>> ensembl = new ArrayList<List<Integer>>();
		List<Integer> mapCPs = new ArrayList<Integer>();
		List<Double> mapMLEs = new ArrayList<Double>();
		List<Double> mapLs = new ArrayList<Double>();
		int first = 0;
		int j=0;
		int last = map_j[j][0].getMaxObject();
		//		mapMLEs.add(segment_mle[first][last-1]);
		j++;

		while(j<m+1){
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getMaxObject()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getMaxObject()-1]));
			//			System.out.println(Util.list(map[last].getMaxObject(),n,map[last].getMaxObjectValue()));
			mapCPs.add(last-1);
			first = last;
			last = map_j[j][last].getMaxObject();
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
		//		System.out.println(Util.list(sample_n(new double[]{.2,.3,.1,.049,.001,.35}, 1000)));
		//		System.exit(0);
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

		int l=40;
		double[][] Y2 = new double[2][l*3];
		double[] Y = new double[l*4];
		MapList<String, Double[]> data = new MapList<String, Double[]>();
		//		RandomEngine rng =  new DRand();
		RandomEngine rng = RandomEngine.makeDefault();

		List<Double> Ys = new ArrayList<Double>();


		int R=250;
		double bbb = .3;
		double ddd = .1;
		NegativeBinomial pd = new NegativeBinomial(R, bbb, rng);

		//		NormalDistributionImpl nd = new NormalDistributionImpl(0, 1);
		for(int i=0; i<l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, ddd, rng);
		for(int i=l; i<2*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, ddd, rng);
		for(int i=2*l; i<3*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, bbb, rng);

		//		NormalDistributionImpl nd = new NormalDistributionImpl(0, 1);
		for(int i=0; i<l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, ddd, rng);
		for(int i=l; i<2*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(R, ddd, rng);
		for(int i=2*l; i<3*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
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

		{
			JFreeChart jfc0 = ChartUtils.createScatterChart(data.getMap(), "", "", "", false);
			ChartUtils.setDomainAxisRange(jfc0, 0, Ys.size());
			ChartUtils.showChart(jfc0);
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
		SegmentationResult sr = doubly_constrained_segmentation(Y, alpha_0, beta_0, nb_r, R, .99, constrained_decreasing, .9);
		JointSegmentationResult jsr = doubly_constrained_multi_segmentation(Y2, alpha_0, beta_0, nb_r, 10, .99, constrained_decreasing, .9);
		PrefixSuffixResult prefix_suffix = prefix_suffix(Y2, alpha_0, beta_0, R, 10, .99);
		
		ConfidenceCalculator cc = new ConfidenceCalculator(prefix_suffix);
		
		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-1, l).log_confidence), cc.calculateConfidence(l-1, l).log_odds);
		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-1, l+1).log_confidence), cc.calculateConfidence(l-1, l+1).log_odds);
		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-5, l+5).log_confidence), cc.calculateConfidence(l-5, l+5).log_odds);
		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-10, l+10).log_confidence), cc.calculateConfidence(l-10, l+10).log_odds);
		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-20, l+20).log_confidence), cc.calculateConfidence(l-20, l+20).log_odds);
		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-40, l+40).log_confidence), cc.calculateConfidence(l-40, l+40).log_odds);
//		System.out.printf("cc %.6f %.3f\n", Math.exp(cc.calculateConfidence(l-80, l+80).log_confidence), cc.calculateConfidence(l-80, l+80).log_odds);
		
		{
			MapList<String, Double[]> ml = new MapList<String, Double[]>();
			for(int i=0;i<jsr.change_point.length;i++){
				ml.put("", new Double[]{jsr.change_point[i]+0.,0.});
				ml.put("", new Double[]{jsr.change_point[i]+0.,1.});
				ml.put("", new Double[]{jsr.change_point[i]+0.,null});
			}
			JFreeChart jfc = ChartUtils.createLineChart(ml.getMap(), "", "", "", false);
			ChartUtils.showChart(jfc);
		}


	}
}
