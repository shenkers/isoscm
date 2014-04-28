package tools;

import java.util.Map;

public class PositionClustering {
	
	public static class Clusterer {
		public GenomicIntervalSet clusters;
		int usW;
		int dsW;
		
		public Clusterer(int width) {
			usW = width/2;
			dsW = width-usW;
			clusters = new GenomicIntervalSet();
		}
		
		public void addPosition(String chr, int position){
			clusters.add(chr, position-usW, position+dsW-1);
		}
	}
	
	public static class StrandedClusterer {
		public StrandedGenomicIntervalSet clusters;
		int usW;
		int dsW;
		
		public StrandedClusterer(int width) {
			usW = width/2;
			dsW = width-usW;
			clusters = new StrandedGenomicIntervalSet();
		}
		
		public void addPosition(String chr, int position, char strand){
			clusters.add(chr, position-usW, position+dsW-1, strand);
		}
	}
	
	public static GenomicIntervalSet clusteredPositionIntervals(GenomicIntervalTree<Map<String,Object>> positions, int width){
		int usW = width/2;
		int dsW = width-usW;
		GenomicIntervalSet clusters = new GenomicIntervalSet();
		for(AnnotatedRegion r : positions){
			clusters.add(r.chr, r.start-usW, r.start+dsW-1);
		}
		
		return clusters;
	}
	
	public static StrandedGenomicIntervalSet clusteredPositionIntervals(StrandedGenomicIntervalTree<Map<String,Object>> positions, int width){
		int usW = width/2;
		int dsW = width-usW;
		StrandedGenomicIntervalSet clusters = new StrandedGenomicIntervalSet();
		for(AnnotatedRegion r : positions){
			clusters.add(r.chr, r.start-usW, r.start+dsW-1,r.strand);
		}
		
		return clusters;
	}
	
	public static void main(String[] args) {
		GenomicIntervalTree<Map<String,Object>> positions = new GenomicIntervalTree<Map<String,Object>>();
		
		positions.add("1", 1, 1);
		positions.add("1", 1, 1);
		positions.add("1", 3, 3);
		positions.add("1", 7, 7);
		positions.add("1", 9, 9);
		
		StrandedGenomicIntervalTree<Map<String,Object>> Spositions = new StrandedGenomicIntervalTree<Map<String,Object>>();
		
		Spositions.add("1", 1, 1,'+');
		Spositions.add("1", 1, 1,'-');
		Spositions.add("1", 3, 3,'+');
		Spositions.add("1", 7, 7,'+');
		Spositions.add("1", 9, 9,'+');
		
		for(AnnotatedRegion r : clusteredPositionIntervals(positions, 4)){
			System.out.println(r);
		}
		
		for(AnnotatedRegion r : clusteredPositionIntervals(Spositions, 2)){
			System.out.println(r+" "+r.strand);
		}
	}

}
