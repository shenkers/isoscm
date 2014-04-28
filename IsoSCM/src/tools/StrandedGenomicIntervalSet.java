package tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import util.Util;

public class StrandedGenomicIntervalSet implements Iterable<AnnotatedRegion>, Serializable{

	private static final long serialVersionUID = 1L;
	//	Map<String, IntervalSet> genomic_intervals;
	Map<Character,GenomicIntervalSet> genomic_intervals;
	
	public StrandedGenomicIntervalSet() {
		genomic_intervals = new HashMap<Character, GenomicIntervalSet>();
	}

	public void remove(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(strand))
			genomic_intervals.put(strand, new GenomicIntervalSet());
		genomic_intervals.get(strand).remove(chr, start, end);
	}

	public void add(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(strand))
			genomic_intervals.put(strand, new GenomicIntervalSet());
		genomic_intervals.get(strand).add(chr, start, end);
	}

	public Iterator<AnnotatedRegion> iterator() {
		List<StrandedGenomicIntervalIteratable> iterables = new LinkedList<StrandedGenomicIntervalSet.StrandedGenomicIntervalIteratable>();
		for(Character strand : genomic_intervals.keySet()){
			iterables.add(new StrandedGenomicIntervalIteratable(genomic_intervals.get(strand), strand));
		}
		return new Util.SequentialIterator<AnnotatedRegion>(iterables).iterator();
	}
	
	public Iterable<AnnotatedRegion> overlappingRegions(final String chr, final int start, final int end, char strand){
		if(genomic_intervals.containsKey(strand))
			return genomic_intervals.get(strand).overlappingRegions(chr, start, end, strand);
		else
			return new ArrayList<AnnotatedRegion>(0);
	}

	class StrandedGenomicIntervalIteratable implements Iterable<AnnotatedRegion>{

		GenomicIntervalSet gis; 
		char strand;

		public StrandedGenomicIntervalIteratable(GenomicIntervalSet gis, char strand) {
			this.gis = gis;
			this.strand = strand;
		}
		

		public Iterator<AnnotatedRegion> iterator() {
			return new StrandedGenomicIntervalIterator(gis,strand);
		}

		class StrandedGenomicIntervalIterator implements Iterator<AnnotatedRegion>{

			Iterator<AnnotatedRegion> gisIterator; 
			char strand;

			public StrandedGenomicIntervalIterator(GenomicIntervalSet gis, char strand) {
				gisIterator = gis.iterator();
				this.strand = strand;
			}

			public boolean hasNext() {
				return gisIterator.hasNext();
			}

			public AnnotatedRegion next() {
				AnnotatedRegion r =  gisIterator.next();
				r.strand = strand;
				return r;
			}

			public void remove() {

			}
		}


	}
	
	public int size() {
		int size = 0;
		for(Character c : genomic_intervals.keySet()){
			size += genomic_intervals.get(c).size();
		}
		return size;
	}


	public static void main(String[] args) {
		StrandedGenomicIntervalSet gis = new StrandedGenomicIntervalSet();
		gis.add("chr1", 21, 23, '-');
		gis.add("chr1", 1, 10, '+');
		gis.add("chr2", 1, 10, '+');
		gis.add("chr1", 20, 30, '+');
		gis.add("chr1", 18, 20, '-');
		
		System.out.println();

		for(AnnotatedRegion ar : gis){
			System.out.println(ar+" "+ar.strand);
		}
	}

	
}
