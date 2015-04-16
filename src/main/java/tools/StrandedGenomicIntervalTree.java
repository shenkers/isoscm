package tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import util.Util;

public class StrandedGenomicIntervalTree<S> implements Iterable<AnnotatedRegion>, Serializable{

	private static final long serialVersionUID = 1L;

	Map<Character, GenomicIntervalTree<S>> genomic_intervals;

	public StrandedGenomicIntervalTree() {
		genomic_intervals = new HashMap<Character, GenomicIntervalTree<S>>();
	}

	/**
	 * assumes that r.start == r.end
	 * @param r
	 */
	public void add(AnnotatedRegion r){
		if(!genomic_intervals.containsKey(r.strand))
			genomic_intervals.put(r.strand, new GenomicIntervalTree<S>());
		GenomicIntervalTree<S> is = genomic_intervals.get(r.strand);
		is.add(r.chr, r.start, r.end);
	}

	public void add(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(strand))
			genomic_intervals.put(strand, new GenomicIntervalTree<S>());
		GenomicIntervalTree<S> is = genomic_intervals.get(strand);
		is.add(chr, start, end);
	}

	public void add(String chr, int start, int end, char strand, S data){
		if(!genomic_intervals.containsKey(strand))
			genomic_intervals.put(strand, new GenomicIntervalTree<S>());
		GenomicIntervalTree<S> is = genomic_intervals.get(strand);
		is.add(chr, start, end, data);
	}

	public void remove(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(strand))
			genomic_intervals.put(strand, new GenomicIntervalTree<S>());
		GenomicIntervalTree<S> is = genomic_intervals.get(strand);
		is.remove(chr,start, end);
	}

	//	public AnnotatedRegion getClosest(String chr, int position, char strand){
	//		if(!genomic_positions.containsKey(strand)){
	//			return null;
	//		}
	//		
	//		GenomicPositionTree is = genomic_positions.get(strand);
	//		
	//		return is.getClosest(chr, position);
	//	}
	//	
	//	/**
	//	 * 
	//	 * @param chr
	//	 * @param position
	//	 * @param strand
	//	 * @return closest upstream element or null if none exists
	//	 */
	//	public AnnotatedRegion getClosestUpstream(String chr, int position, char strand){
	//		if(!genomic_positions.containsKey(strand)){
	//			return null;
	//		}
	//		
	//		GenomicPositionTree is = genomic_positions.get(strand);
	//		
	//		if(strand == '+'){
	//			AnnotatedRegion ar = is.getClosestUpstream(chr, position);
	//			if(ar!=null)
	//			ar.strand=strand;
	//			return ar;
	//		}
	//		else{
	//			AnnotatedRegion ar = is.getClosestDownstream(chr, position);
	//			if(ar!=null)
	//				ar.strand=strand;
	//			return ar;
	//		}
	//		}
	//	
	//	/**
	//	 * 
	//	 * @param chr
	//	 * @param position
	//	 * @param strand
	//	 * @return closest downstream or null if none exists
	//	 */
	//	public AnnotatedRegion getClosestDownstream(String chr, int position, char strand){
	//		if(!genomic_positions.containsKey(strand)){
	//			return null;
	//		}
	//		
	//		GenomicPositionTree is = genomic_positions.get(strand);
	//		
	//		if(strand == '+')
	//			return is.getClosestDownstream(chr, position);
	//		else
	//			return is.getClosestUpstream(chr, position);
	//	}
	//
	public boolean contains(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(strand))
			return false;
		else{
			GenomicIntervalTree<S> is = genomic_intervals.get(strand);
			return is.contains(chr, start, end);
		}
	}
	
	public AnnotatedRegion get(String chr, int start, int end, char strand){
		if(!genomic_intervals.containsKey(strand))
			return null;
		else{
			GenomicIntervalTree<S> is = genomic_intervals.get(strand);
			return is.get(chr, start, end, strand);
		}
	}

	public AnnotatedRegion getClosest(String chr, int position, char strand){
		if(!genomic_intervals.containsKey(strand)){
			return null;
		}

		GenomicIntervalTree<S> is = genomic_intervals.get(strand);

		AnnotatedRegion ar = is.getClosest(chr, position);
		
		if(ar!=null)
			ar.strand=strand;

		return ar;
	}

	/**
	 * 
	 * @param chr
	 * @param position
	 * @param strand
	 * @return closest upstream element with start less than or equal to 'position' or null if none exists
	 */
	public AnnotatedRegion getClosestUpstream(String chr, int position, char strand){
		if(!genomic_intervals.containsKey(strand)){
			return null;
		}

		GenomicIntervalTree<S> is = genomic_intervals.get(strand);
		
		AnnotatedRegion ar = null;

		if(strand == '+'){
			ar = is.getClosestUpstream(chr, position);
		}
		else if(strand == '-'){
			ar = is.getClosestDownstream(chr, position);
		
		}
		
		if(ar!=null)
			ar.strand=strand;

		return ar;
	}

	/**
	 * 
	 * @param chr
	 * @param position
	 * @param strand
	 * @return closest downstream or null if none exists
	 */
	public AnnotatedRegion getClosestDownstream(String chr, int position, char strand){
		if(!genomic_intervals.containsKey(strand)){
			return null;
		}

		GenomicIntervalTree<S> is = genomic_intervals.get(strand);

		AnnotatedRegion ar = null;

		if(strand == '+'){
			ar = is.getClosestDownstream(chr, position);
		}else if(strand == '-'){
			ar = is.getClosestUpstream(chr, position);
		}

		if(ar!=null)
			ar.strand=strand;

		return ar;
	}

	public Iterator<AnnotatedRegion> iterator() {
		List<StrandedGenomicIntervalIteratable> iterables = new LinkedList<StrandedGenomicIntervalIteratable>();
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

	//	public Iterable<IntervalNode> overlappingNodes(final String chr, final int start, final int end, char strand){
	//		if(genomic_intervals.containsKey(strand))
	//			return genomic_intervals.get(strand).overlappingNodes(chr, start, end, strand);
	//		else
	//			return new ArrayList<IntervalNode>(0);
	//	}

	class StrandedGenomicIntervalIteratable implements Iterable<AnnotatedRegion>{

		GenomicIntervalTree<S> gis; 
		char strand;

		public StrandedGenomicIntervalIteratable(GenomicIntervalTree<S> genomicIntervalTree, char strand) {
			this.gis = genomicIntervalTree;
			this.strand = strand;
		}


		public Iterator<AnnotatedRegion> iterator() {
			return new StrandedGenomicIntervalIterator(this.gis,this.strand);
		}

		class StrandedGenomicIntervalIterator implements Iterator<AnnotatedRegion>{

			Iterator<AnnotatedRegion> gisIterator; 
			char strand;

			public StrandedGenomicIntervalIterator(GenomicIntervalTree<S> gis, char strand) {
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

	public int size(){
		int n=0;

		for(char strand : genomic_intervals.keySet()){
			n += genomic_intervals.get(strand).size();
		}

		return n;
	}

	public void clear(){
		genomic_intervals.clear();
	}
	
	public static void main(String[] args) {
		StrandedGenomicIntervalTree<Object> gis = new StrandedGenomicIntervalTree();

		gis.add("chr1", 21, 23, '-');
		gis.add("chr1", 1, 10, '+');
		gis.add("chr2", 1, 10, '+');
		gis.add("chr1", 20, 30, '+');
		gis.add("chr1", 18, 20, '-');

		System.out.println();

		for(AnnotatedRegion ar : gis){
			System.out.println(ar+" "+ar.strand);
		}

		System.out.println();

		for(AnnotatedRegion r : gis.overlappingRegions("chr1", 10, 20,'+')){
			System.out.println(r);
		}
		System.out.println("done");

	}

}
