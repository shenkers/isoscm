package scm;

public class PrefixSuffixResult {
	// likelihood of the sequence of observations between indexes t and s
	double[][] log_p;
	// likelihood of the suffix starting at position i
	double[] log_q_suffix;
	// likelihood of the prefix ending at position i
	double[] log_q_prefix;
	// likelihood of the prefix ending at position i, given no change-points
	double[] log_prefix;
	// likelihood of the suffix starting at position i, given no change-points
	double[] log_suffix;

	
	// length probability density function of segment length
	double[] log_g;
	// length cumulative density function of segment length
	double[] log_G;
	// length of the sequence of observations
	int l;
	public PrefixSuffixResult(double[][] log_p, double[] log_g, double[] log_G,
			double[] log_q_prefix, double[] log_q_suffix, double[] log_prefix,
			double[] log_suffix, int l) {
		this.log_p = log_p;
		this.log_g = log_g;
		this.log_G = log_G;
		this.log_q_prefix = log_q_prefix;
		this.log_q_suffix = log_q_suffix;
		this.log_prefix = log_prefix;
		this.log_suffix = log_suffix;
		this.l = l;
	}
	
	
}
