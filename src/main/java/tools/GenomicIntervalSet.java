package tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class GenomicIntervalSet implements Iterable<AnnotatedRegion>, Serializable{

	private static final long serialVersionUID = 1L;

	Map<String, IntegerIntervalSet> genomic_intervals;

	public GenomicIntervalSet() {
		genomic_intervals = new HashMap<String, IntegerIntervalSet>();
	}

	public void remove(String chr, int start, int end){
		if(!genomic_intervals.containsKey(chr))
			genomic_intervals.put(chr, new IntegerIntervalSet());
		IntegerIntervalSet is = genomic_intervals.get(chr);
		is.remove(start, end);
	}

	public void add(String chr, int start, int end){
		if(!genomic_intervals.containsKey(chr))
			genomic_intervals.put(chr, new IntegerIntervalSet());
		IntegerIntervalSet is = genomic_intervals.get(chr);
		is.add(start, end);
	}

	public Iterator<AnnotatedRegion> iterator() {
		return new GenomicIntervalIterator(genomic_intervals);
	}

	public int size(){
		int n=0;

		for(String chr : genomic_intervals.keySet()){
			n += genomic_intervals.get(chr).size();
		}

		return n;
	}

	public Iterable<AnnotatedRegion> overlappingRegions(final String chr, final int start, final int end){
		return new OverlappingRegionsIterable(chr, start, end, '.');
	}
	
	public Iterable<AnnotatedRegion> overlappingRegions(final String chr, final int start, final int end, char strand){
		return new OverlappingRegionsIterable(chr, start, end, strand);
	}

	class OverlappingRegionsIterable implements Iterable<AnnotatedRegion>{

		String chr;
		int start,end;
		char strand;

		public OverlappingRegionsIterable(String chr, int start, int end, char strand) {
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.strand = strand;
		}

		public Iterator<AnnotatedRegion> iterator() {
			if(genomic_intervals.containsKey(chr)){
				final IntegerIntervalSet is = genomic_intervals.get(chr);
				final Iterator<List<Integer>> it = is.overlappingIntervalIterable(start, end).iterator();
				Iterator<AnnotatedRegion> iar = new Iterator<AnnotatedRegion>() {
					public void remove() {
					}
					public AnnotatedRegion next() {
						List<Integer> next = it.next();
						return new AnnotatedRegion("", chr, next.get(0), next.get(1), strand);
					}
					public boolean hasNext() {
						return it.hasNext();
					}
				};
				return iar;
			}
			else{
				return new ArrayList<AnnotatedRegion>(0).iterator();
			}
		}
	}

	class GenomicIntervalIterator implements Iterator<AnnotatedRegion>{

		Map<String, IntegerIntervalSet> genomic_intervals;
		Iterator<String> chr_iterator;
		IntegerIntervalSet is;
		boolean has_next;
		IntervalNode node;
		String chr;

		public GenomicIntervalIterator(Map<String, IntegerIntervalSet> genomic_intervals) {
			this.genomic_intervals = genomic_intervals;
			chr_iterator = genomic_intervals.keySet().iterator();
		
			if(chr_iterator.hasNext()){
				increment_chr();
			}
			else{
				has_next = false;
			}
		}

		private void increment_chr(){
			has_next=false;
			while(chr_iterator.hasNext()){
				chr = chr_iterator.next();
				is = genomic_intervals.get(chr);
				node = is.tree.getSmallest();
				has_next = node!=null;
				if(has_next)
					break;
			}
		}

		public boolean hasNext() {
			return has_next;
		}

		public AnnotatedRegion next() {
			AnnotatedRegion r = new AnnotatedRegion("", chr, node.start, node.end, '.');
			node = node.getNext();
			if(node==null){
				increment_chr();
			}
			return r;
		}

		public void remove() {

		}

	}

	public static void main(String[] args) {
		GenomicIntervalSet gis = new GenomicIntervalSet();
		gis.add("chr1", 1, 10);
		gis.add("chr1", -12, 16);
		gis.add("chr2", 1, 10);
		gis.add("chr1", 20, 30);
		gis.add("chr1", 18, 19);

		for(AnnotatedRegion ar : gis){
			System.out.println(ar);
		}

		System.out.println(gis.size());
	}

}
