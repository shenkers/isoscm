package scm;

public class ConfidenceResult {

	double log_confidence;
	double log_odds;
	
	public ConfidenceResult(double confidence, double log_odds) {
		this.log_confidence = confidence;
		this.log_odds = log_odds;
	}	
}
