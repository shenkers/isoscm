package scm;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;

import util.Util;
import util.Util.ExtremeObjectTracker;
import util.Util.IntegralTable;
import util.Util.MapList;
import cern.jet.random.NegativeBinomial;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;

public class BetaNegativeBinomial {

	

	/*
	 * Calculate the probability of every sequence from (s,t) of observations
	 * y coming from a distribution with a single parameter, integrating
	 * over possible parameter values with prior alpha and beta. 
	 */
	public static SCMResult calculate_segment_probabilities(double[] y, double alpha_0, double beta_0, int r){	
	
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
	public static SegmentationResult viterbi_segmentation(double[] y, double alpha_0, double beta_0, int nb_r, int r, double p){
		double[] log_q	= new double[y.length];
		@SuppressWarnings("unchecked")
		ExtremeObjectTracker<Integer, Double>[] map = new ExtremeObjectTracker[y.length];

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
				double p_segment_length = t==0? Math.log(1 - Math.exp(log_G[s])) : log_g[s-t];
				map[t].put(s+1, log_p[t][s]+map[s+1].getExtremeValue()+ p_segment_length);			

				// sum over the positions of the next change point
				log_probs.add(log_p[t][s]+log_q[s+1]+log_g[s-t]);
			}

			// the probability of there being no more changepoints
			map[t].put(n, log_p[t][n-1] + Math.log(1 - Math.exp(log_G[n-t-1])));
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
			//			System.out.println(Util.list(Util.list(first,last-1),Util.list(last,map[last].getExtreme()-1)));
			//			System.out.println(Util.list(segment_mle[first][last-1],segment_mle[last][map[last].getExtreme()-1]));
			//			System.out.println(Util.list(map[last].getExtreme(),n,map[last].getExtremeValue()));
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
//		JFreeChart jfc = ChartUtils.createLineChart(dat.getMap(), "", "", "", false);
//		ChartUtils.showChart(jfc);

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

	public static void main(String[] args) throws MathException {
	
//		testFlyData();
//				testRealData();

		int l=50;
		double[][] Y2 = new double[2][l*3];
		double[] Y = new double[l*3];
		MapList<String, Double[]> data = new MapList<String, Double[]>();
		RandomEngine rng = new DRand();
		int nb_r = 100;
		NegativeBinomial pd = new NegativeBinomial(nb_r, .25, rng);
//		NormalDistributionImpl nd = new NormalDistributionImpl(0, 1);
		for(int i=0; i<l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(nb_r, .01, rng);
		for(int i=l; i<2*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		pd = new NegativeBinomial(nb_r, .25, rng);
		for(int i=2*l; i<3*l; i++){
			Y[i] = pd.nextDouble();
			Y2[0][i] = Y[i];
			data.put("1", new Double[]{(double) i,Y[i]});
		}
		l=(3*l)/2;
		pd = new NegativeBinomial(100, .15, rng);
		for(int i=0; i<1*l; i++){
			Y2[1][i] = pd.nextDouble();
//			data.put("2", new Double[]{(double) i,Y2[1][i]});
		}
		pd = new NegativeBinomial(100, .01, rng);
		for(int i=1*l; i<2*l; i++){
			Y2[1][i] = pd.nextDouble();
//			data.put("2", new Double[]{(double) i,Y2[1][i]});
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
		
//		SegmentationResult sr = viterbi_segmentation(Y2, alpha_0, beta_0, nb_r, 1, .95, 2);
		SegmentationResult sr = viterbi_segmentation(Y, alpha_0, beta_0, nb_r, 2, .99);
	
//		SegmentationResult sr = viterbi_segmentation(Y, m_0, V_0, a_0, b_0, 1, .95);
//		System.out.println(Util.list(sr.change_point));
		for(int i=0; i<sr.change_point.length; i++){
//			data.put("a", new Double[]{sr.change_point[i]+0.,0.});
			data.put("a", new Double[]{sr.change_point[i]+0.,-1.});
		}
		
//		ChartUtils.showChart(ChartUtils.createScatterChart(data.getMap(), "", "", "", false));

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
