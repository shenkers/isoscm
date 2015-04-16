package tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class GenomicIntervalTree<S> implements Iterable<AnnotatedRegion>, Serializable{

	private static final long serialVersionUID = 1L;

	Map<String, IntegerIntervalTree<S>> genomic_intervals;

	public GenomicIntervalTree() {
		genomic_intervals = new HashMap<String, IntegerIntervalTree<S>>();
	}

	public void add(String chr, int start, int end){
		if(!genomic_intervals.containsKey(chr))
			genomic_intervals.put(chr, new IntegerIntervalTree<S>());
		IntegerIntervalTree<S> is = genomic_intervals.get(chr);
		is.insert(start,end);
	}

	public void add(String chr, int start, int end, S data){
		if(!genomic_intervals.containsKey(chr))
			genomic_intervals.put(chr, new IntegerIntervalTree<S>());
		IntegerIntervalTree<S> is = genomic_intervals.get(chr);
		is.insert(start,end,data);
	}

	public void remove(String chr, int start, int end){
		if(!genomic_intervals.containsKey(chr))
			genomic_intervals.put(chr, new IntegerIntervalTree<S>());
		IntegerIntervalTree<S> is = genomic_intervals.get(chr);
		is.delete(start, end);
	}

	public boolean contains(String chr, int start, int end){
		if(!genomic_intervals.containsKey(chr))
			return false;
		else{
			IntegerIntervalTree<S> is = genomic_intervals.get(chr);
			return is.contains(start, end);
		}
	}
	
	public AnnotatedRegion get(String chr, int start, int end){
		return get(chr, start, end, '.');
	}

	public AnnotatedRegion get(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(chr))
			return null;
		else{
			IntegerIntervalTree<S> is = genomic_intervals.get(chr);
			IntervalNode<S> n = is.get(start, end);
			AnnotatedRegion r = new AnnotatedRegion("", chr, n.getStart(), n.getEnd(), strand);
			r.attributes = (Map) n.data;
			return r;
		}
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

	//	public Iterable<IntervalNode> overlappingNodes(final String chr, final int start, final int end){
	//		return new OverlappingNodeIterable(chr, start, end, '.');
	//	}

	public Iterable<AnnotatedRegion> overlappingRegions(final String chr, final int start, final int end, char strand){
		return new OverlappingRegionsIterable(chr, start, end, strand);
	}

	//	public Iterable<IntervalNode> overlappingNodes(final String chr, final int start, final int end, char strand){
	//		return new OverlappingNodeIterable(chr, start, end, strand);
	//	}
	//	
	//	class OverlappingNodeIterable implements Iterable<IntervalNode>{
	//
	//		String chr;
	//		int start,end;
	//		char strand;
	//
	//		public OverlappingNodeIterable(String chr, int start, int end, char strand) {
	//			this.chr = chr;
	//			this.start = start;
	//			this.end = end;
	//			this.strand = strand;
	//		}
	//
	//		public Iterator<IntervalNode> iterator() {
	//			if(genomic_intervals.containsKey(chr)){
	//				final IntegerIntervalTree is = genomic_intervals.get(chr);
	//				final Iterator<IntervalNode> it = is.overlaps(start, end).iterator();
	//				Iterator<IntervalNode> iar = new Iterator<IntervalNode>() {
	//					public void remove() {
	//					}
	//					public IntervalNode next() {
	//						IntervalNode next = it.next();
	//						return next;
	//					}
	//					public boolean hasNext() {
	//						return it.hasNext();
	//					}
	//				};
	//				return iar;
	//			}
	//			else{
	//				return new ArrayList<IntervalNode>(0).iterator();
	//			}
	//		}
	//	}

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
				final IntegerIntervalTree<S> is = genomic_intervals.get(chr);
				final Iterator<IntervalNode<S>> it = is.overlaps(start, end).iterator();
				Iterator<AnnotatedRegion> iar = new Iterator<AnnotatedRegion>() {
					public void remove() {
					}
					public AnnotatedRegion next() {
						IntervalNode<S> next = it.next();
						AnnotatedRegion r = new AnnotatedRegion("", chr, next.getStart(), next.getEnd(), strand);
						r.attributes = (Map) next.data;
						return r;
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

		Map<String, IntegerIntervalTree<S>> genomic_intervals;
		Iterator<String> chr_iterator;
		IntegerIntervalTree<S> is;
		IntervalNode<S> in;
		boolean has_next;
		String chr;

		public GenomicIntervalIterator(Map<String, IntegerIntervalTree<S>> genomic_intervals) {
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
				in = is.getSmallest();
				has_next = in!=null;
				if(has_next)
					break;
			}
		}

		public boolean hasNext() {
			return has_next;
		}

		public AnnotatedRegion next() {
			AnnotatedRegion r = new AnnotatedRegion("", chr, in.getStart(), in.getEnd(), '.', (Map<String, Object>) in.data);
			in = in.getNext();

			if(in==null){
				increment_chr();
			}
			return r;
		}

		public void remove() {

		}

	}

	public AnnotatedRegion getClosest(String chr, int position){
		if(!genomic_intervals.containsKey(chr))
			return null;

		IntegerIntervalTree<S> is = genomic_intervals.get(chr);

		IntervalNode<S> n = is.getNotLarger(position);

		if(n==null)
			n = is.getNotSmaller(position);
		if(n==null)
			return null;

		int closestPosition = n.start;
		S data = n.data;
		n = n.getNext();
		if(n!=null){
		int d_position = closestPosition-position;
		int d_next = n.start-position;
		while(n!=null && Math.abs(d_next) <= Math.abs(d_position)){
			closestPosition = n.start;
			d_position = closestPosition-position;
			data = n.data;
			n = n.getNext();
			if(n!=null)
			d_next = n.start-position;
		}
		}

		return new AnnotatedRegion("", chr, closestPosition, closestPosition, '.', (Map) data);
	}


	public AnnotatedRegion getClosestUpstream(String chr, int position){
		if(!genomic_intervals.containsKey(chr))
			return null;

		IntegerIntervalTree<S> is = genomic_intervals.get(chr);

		IntervalNode<S> n = is.getNotLarger(position);

		if(n==null)
			return null;

		int closestPosition = n.start;
		S data = n.data;

		return new AnnotatedRegion("", chr, closestPosition, closestPosition, '.', (Map) data);
	}

	public AnnotatedRegion getClosestDownstream(String chr, int position){
		if(!genomic_intervals.containsKey(chr))
			return null;

		IntegerIntervalTree<S> is = genomic_intervals.get(chr);

		IntervalNode<S> n = is.getNotSmaller(position);

		if(n==null)
			return null;

		int closestPosition = n.start;
		S data = n.data;

		return new AnnotatedRegion("", chr, closestPosition, closestPosition, '.', (Map) data);
	}

	public static void main(String[] args) {
		GenomicIntervalTree<Object> gis = new GenomicIntervalTree();

		gis.add("chr1", 1, 10);
		gis.add("chr1", 1, 10);
		gis.add("chr1", 5, 10);
		gis.add("chr1", 15, 20);
		gis.add("chr2", 15, 20);

		for(AnnotatedRegion r : gis){
			System.out.println(r);
		}
		System.out.println();
		for(AnnotatedRegion r : gis.overlappingRegions("chr2", 20, 25)){
			System.out.println(r);
		}

		System.out.println();
		System.out.println(gis.getClosestDownstream("chr1", -1));
		System.out.println(gis.getClosestDownstream("chr1", 0));
		System.out.println(gis.getClosestDownstream("chr1", 1));
		System.out.println(gis.getClosestDownstream("chr1", 2));
		System.out.println(gis.getClosestDownstream("chr1", 3));
		System.out.println(gis.getClosestDownstream("chr1", 4));
		System.out.println(gis.getClosestDownstream("chr1", 5));
		System.out.println(gis.getClosestDownstream("chr1", 6));
		System.out.println(gis.getClosestDownstream("chr1", 14));
		System.out.println(gis.getClosestDownstream("chr1", 15));
		System.out.println(gis.getClosestDownstream("chr1", 16));
		}

	public void clear() {
		genomic_intervals.clear();
	}

}
