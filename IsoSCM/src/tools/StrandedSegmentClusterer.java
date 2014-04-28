package tools;

import java.util.Iterator;

public class StrandedSegmentClusterer implements Iterable<AnnotatedRegion>{

	public StrandedGenomicIntervalSet clusters;
	int usW;
	int dsW;

	public StrandedSegmentClusterer(int width) {
		usW = width/2;
		dsW = Math.max(width-usW,0);
		clusters = new StrandedGenomicIntervalSet();
	}

	public void addSegment(String chr, int start, int end, char strand){
		clusters.add(chr, start-usW, end+dsW, strand);
	}
	
	public void addSegment(String chr, int start, int end, char strand, boolean extra_width){
		if(extra_width)
			addSegment(chr, start, end, strand);
		else
			clusters.add(chr, start, end, strand);
	}

	class MergedSegmentIterator implements Iterator<AnnotatedRegion>{

		Iterator<AnnotatedRegion> iterator;

		public MergedSegmentIterator(Iterator<AnnotatedRegion> iterator) {
			this.iterator = iterator;
		}

		public boolean hasNext() {
			return iterator.hasNext();
		}

		public AnnotatedRegion next() {
			AnnotatedRegion r =  iterator.next();
			r.start = r.start + usW;
			r.end = r.end - dsW;
			return r;
		}

		public void remove() {

		}
	}
	
	public Iterator<AnnotatedRegion> iterator() {
		return new MergedSegmentIterator(clusters.iterator());
	}
	
	public static void main(String[] args) {
		StrandedSegmentClusterer ssc = new StrandedSegmentClusterer(1);
		ssc.addSegment("", 1, 1, '.');
		ssc.addSegment("", 5, 5, '.');
		ssc.addSegment("", 3, 3, '.');
		
		for(AnnotatedRegion r : ssc){
			System.out.println(r);
		}
	}
}
