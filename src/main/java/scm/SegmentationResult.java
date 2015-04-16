package scm;

public class SegmentationResult {
	// the positions of the change points, each changepoint marks the position
	// of the last observation from a segment
	public int[] change_point;
	// the MLE of each segment, this is equal to the number of change points + 1
	public double[] segment_mle;
	// the log-likelihood of each changepoint
	public double[] change_point_ll;
	// the change points occur on an interval from [0-l) -- including zero, excluding 'l'
	public int l;
	// the number of segments
	public int n_segments;

	public SegmentationResult(int[] change_point, double[] segment_mle, double[] change_point_ll, int l) {
		this.change_point = change_point;
		this.segment_mle = segment_mle;
		this.change_point_ll = change_point_ll;
		this.l = l;
		n_segments = segment_mle.length;			
	}

	public int[][] segments(){
		int[][] segments = new int[n_segments][2];

		if(n_segments==1){
			segments[0] = new int[]{0,l-1};
		}
		else{
			segments[0] = new int[]{0,change_point[0]};
			for(int i=0; i<change_point.length-1; i++){
				segments[i+1] = new int[]{change_point[i]+1,change_point[i+1]};
			}
			segments[change_point.length] = new int[]{change_point[change_point.length-1]+1,l-1};
		}

		return segments;
	}
}
