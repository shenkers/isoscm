package tools;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import util.Util;
import util.Util.Criteria;
import util.Util.ExtremeTracker;
import util.Util.Instantiator;
import util.Util.MapFactory;

public class IntervalTools {

	/**
	 * 
	 * @param ref_position
	 * @param position
	 * @param isNegativeStrand
	 * @param distanceUpstream
	 * @return the distance up(down)stream position is from ref_position, assuming position is on the given strand
	 */
	public static int offsetAmount(int ref_position, int position, boolean isNegativeStrand, boolean distanceUpstream){
		int i = isNegativeStrand ^ distanceUpstream ? 1 : -1;
		return (ref_position-position)*i;
	}

	public static int[] offsetAmounts(int ref_start, int ref_end, int start_position, int end_position){
		return new int[]{(ref_start-start_position),-(ref_end-end_position)};
	}

	public static int offsetPosition(int position, int amount, boolean isNegativeStrand, boolean upstream){
		int i = isNegativeStrand ^ upstream ? -1 : 1;
		return position + (amount*i);
	}

	public static int[] offsetInterval(int position, int upstream_amount, int downstream_amount, boolean isNegativeStrand){
		int s = isNegativeStrand ? downstream_amount :  upstream_amount;
		int e = !isNegativeStrand ? downstream_amount :  upstream_amount;
		return new int[]{position - s, position + e};
	}

	public static int[] offsetInterval(int start, int end, int upstream_amount, int downstream_amount, boolean isNegativeStrand){
		int s = isNegativeStrand ? downstream_amount :  upstream_amount;
		int e = !isNegativeStrand ? downstream_amount :  upstream_amount;
		return new int[]{start - s, end + e};
	}

	/**
	 * 
	 * @param s
	 * @param e
	 * @param start
	 * @param eend
	 * @return whether (s,e) is nested within (start,end)
	 */
	public static boolean isContained(int s, int e, int start, int end){
		return s> start-1 && e < end+1 && s < e+1;
	}

	/**
	 * 
	 * @param s
	 * @param e
	 * @param start
	 * @param eend
	 * @return whether (s,e) is nested within (start,end)
	 */
	public static boolean overlaps(int s, int e, int start, int end){
		return (s> start-1 && s < end+1) || (e> start-1 && e < end+1) || (start> s-1 && end < e+1);
	}

	public static int intersectionLength(int s1, int e1, int s2, int e2){
		return Math.min(e1, e2) - Math.max(s1,s2)+1;
	}

	public static int[] intersection(int s1, int e1, int s2, int e2){
		if(!overlaps(s1, e1, s2, e2))
			throw new RuntimeException(Util.sprintf("intervals (%d,%d) and (%d,%d) do not overlap",s1,e1,s2,e2));
		return new int[]{Math.max(s1,s2), Math.min(e1, e2)};
	}

	public static int[] union(int s1, int e1, int s2, int e2){
		if(!overlaps(s1-1, e1+1, s2, e2))
			throw new RuntimeException(Util.sprintf("intervals (%d,%d) and (%d,%d) do not overlap",s1,e1,s2,e2));
		return new int[]{Math.min(s1,s2), Math.max(e1, e2)};
	}

	public static class OverlappingIterable implements Iterable<AnnotatedRegion>{

		Iterable<AnnotatedRegion> iterable;
		StrandedGenomicIntervalTree<Map<String,Object>> regions;

		public OverlappingIterable(Iterable<AnnotatedRegion> iterable, StrandedGenomicIntervalTree<Map<String,Object>> regions){
			this.iterable=iterable;
			this.regions=regions;
		}

		public Iterator<AnnotatedRegion> iterator() {
			return new OverlappingIterator(iterable.iterator(), regions);
		}

	}

	static class OverlappingIterator implements Iterator<AnnotatedRegion>{

		Iterator<AnnotatedRegion> it;
		AnnotatedRegion next;
		boolean hasNext;
		StrandedGenomicIntervalTree<Map<String,Object>> regions;

		public OverlappingIterator(	Iterator<AnnotatedRegion> it, StrandedGenomicIntervalTree<Map<String,Object>> regions) {
			this.it = it;
			this.regions=regions;
			findNext();
		}

		public boolean hasNext() {
			return hasNext;
		}

		void findNext(){
			while(it.hasNext()){
				AnnotatedRegion n = it.next();
				if(regions.overlappingRegions(n.chr, n.start, n.end, n.strand).iterator().hasNext()){
					next=n;
					hasNext=true;
					return;
				}
			}
			next=null;
			hasNext=false;
		}

		public AnnotatedRegion next() {
			AnnotatedRegion last = next;
			findNext();
			return last;
		}

		public void remove() {
			// TODO Auto-generated method stub
		}

	}

	public static class AttributeIterable implements Iterable<AnnotatedRegion>{

		Iterable<AnnotatedRegion> iterable;
		String attribute;
		Object value;

		public AttributeIterable(Iterable<AnnotatedRegion> iterable, String attribute, Object value){
			this.iterable=iterable;
			this.attribute=attribute;
			this.value=value;
		}

		public Iterator<AnnotatedRegion> iterator() {
			return new AttributeIterator(iterable.iterator(), attribute, value);
		}

	}

	static class AttributeIterator implements Iterator<AnnotatedRegion>{

		Iterator<AnnotatedRegion> it;
		AnnotatedRegion next;
		boolean hasNext;
		String attribute;
		Object value;

		public AttributeIterator(	Iterator<AnnotatedRegion> it, String attribute, Object value) {
			this.it = it;
			this.attribute=attribute;
			this.value=value;
			findNext();
		}

		public boolean hasNext() {
			return hasNext;
		}

		void findNext(){
			while(it.hasNext()){
				AnnotatedRegion n = it.next();
				if(value.equals(n.getAttribute(attribute))){
					next=n;
					hasNext=true;
					return;
				}
			}
			next=null;
			hasNext=false;
		}

		public AnnotatedRegion next() {
			AnnotatedRegion last = next;
			findNext();
			return last;
		}

		public void remove() {
			// TODO Auto-generated method stub
		}

	}

	public static class AnnotatedRegionCriteriaIterable implements Iterable<AnnotatedRegion>{

		Iterable<AnnotatedRegion> iterable;
		Criteria<AnnotatedRegion> c;

		public AnnotatedRegionCriteriaIterable(Iterable<AnnotatedRegion> iterable, Criteria<AnnotatedRegion> c){
			this.iterable=iterable;
			this.c=c;
		}

		public Iterator<AnnotatedRegion> iterator() {
			return new Util.CriteriaIterator<AnnotatedRegion>(iterable.iterator(), c);
		}

	}



	public static StrandedGenomicIntervalTree<Map<String,Object>> aggregateIntervals(Iterable<AnnotatedRegion> annotations, String groupingAttribute){
		StrandedGenomicIntervalTree<Map<String,Object>> intervals = new StrandedGenomicIntervalTree<Map<String,Object>>();

		MapFactory<Object, ExtremeTracker<Integer>> transcript_termini = new MapFactory<Object, ExtremeTracker<Integer>>(new Instantiator<ExtremeTracker<Integer>>() {
			public ExtremeTracker<Integer> instantiate(Object... objects) {
				return new ExtremeTracker<Integer>();
			}
		});
		Map<Object,Boolean> transcript_isNegativeStrand = new HashMap<Object, Boolean>();
		Map<Object,String> transcript_chr = new HashMap<Object, String>();
		for(AnnotatedRegion r : annotations){
			Object grouping_id = r.getAttribute(groupingAttribute);
			transcript_chr.put(grouping_id, r.chr);
			transcript_isNegativeStrand.put(grouping_id, r.isNegativeStrand());
			ExtremeTracker<Integer> gene_termini = transcript_termini.get(grouping_id);
			gene_termini.put(r.start);
			gene_termini.put(r.end);
		}

		for(Object grouping_id : transcript_termini.getMap().keySet()){
			boolean isNegativeStrand = transcript_isNegativeStrand.get(grouping_id);
			String chr = transcript_chr.get(grouping_id);
			ExtremeTracker<Integer> gene_termini = transcript_termini.get(grouping_id);
			int start = gene_termini.getMin();
			int end = gene_termini.getMax();
			char strand = isNegativeStrand ? '-' : '+';
			Map<String, Object> data = new HashMap<String, Object>();
			data.put(groupingAttribute, grouping_id);
			intervals.add(chr, start, end, strand, data);
		}

		return intervals;
	}

	public static StrandedGenomicIntervalTree<Map<String,Object>> aggregateIntervals(Iterable<AnnotatedRegion> annotations, List<String> attributes, String groupingAttribute){
		StrandedGenomicIntervalTree<Map<String,Object>> intervals = new StrandedGenomicIntervalTree<Map<String,Object>>();

		MapFactory<String, ExtremeTracker<Integer>> transcript_termini = new MapFactory<String, ExtremeTracker<Integer>>(new Instantiator<ExtremeTracker<Integer>>() {
			public ExtremeTracker<Integer> instantiate(Object... objects) {
				return new ExtremeTracker<Integer>();
			}
		});
		Map<String,Boolean> transcript_isNegativeStrand = new HashMap<String, Boolean>();
		Map<String,String> transcript_chr = new HashMap<String, String>();

		MapFactory<String, Map<String,Object>> data = new MapFactory<String, Map<String,Object>>(new Instantiator<Map<String,Object>>() {
			public Map<String, Object> instantiate(Object... objects) {
				return new HashMap<String, Object>();
			}
		});

		for(AnnotatedRegion r : annotations){
			String grouping_id = (String) r.getAttribute(groupingAttribute);
			transcript_chr.put(grouping_id, r.chr);
			transcript_isNegativeStrand.put(grouping_id, r.isNegativeStrand());
			ExtremeTracker<Integer> gene_termini = transcript_termini.get(grouping_id);
			gene_termini.put(r.start);
			gene_termini.put(r.end);
			Map<String, Object> d = data.get(grouping_id);
			for(String attribute : attributes){
				d.put(attribute, r.getAttribute(attribute));
			}
		}

		for(String grouping_id : transcript_termini.getMap().keySet()){
			boolean isNegativeStrand = transcript_isNegativeStrand.get(grouping_id);
			String chr = transcript_chr.get(grouping_id);
			ExtremeTracker<Integer> gene_termini = transcript_termini.get(grouping_id);
			int start = gene_termini.getMin();
			int end = gene_termini.getMax();
			char strand = isNegativeStrand ? '-' : '+';
			intervals.add(chr, start, end, strand, data.get(grouping_id));
		}

		return intervals;
	}

	public static List<AnnotatedRegion> BoundaryIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int pos, char strand, boolean is5p){
		List<AnnotatedRegion> matchingIntervals = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, pos, pos, strand)){
			if((is5p?r.get5Prime():r.get3Prime()) == pos)
				matchingIntervals.add(r);
		}

		return matchingIntervals;
	}

	public static List<AnnotatedRegion> BoundaryIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand, boolean is5p){
		List<AnnotatedRegion> matchingIntervals = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			int p = is5p?r.get5Prime():r.get3Prime();
			if(overlaps(p, p, start, end))
				matchingIntervals.add(r);
		}

		return matchingIntervals;
	}

	/**
	 * 
	 * @param intervals
	 * @param chr
	 * @param start
	 * @param end
	 * @param strand
	 * @return all intervals that match the specified boundary
	 */
	public static List<AnnotatedRegion> MatchingIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> matchingIntervals = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			if(r.start==start && r.end==end)
				matchingIntervals.add(r);
		}

		return matchingIntervals;
	}


	/**
	 * 
	 * @param intervals
	 * @param chr
	 * @param start
	 * @param end
	 * @param strand
	 * @return all intervals have either a 5' or 3' boundary within the specified segment
	 */
	public static List<AnnotatedRegion> BoundaryIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> matchingIntervals = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			if(overlaps(r.start, r.start, start, end)||overlaps(r.end, r.end, start, end))
				matchingIntervals.add(r);
		}

		return matchingIntervals;
	}


	public static List<AnnotatedRegion> SelectIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand, String attribute, Object value){
		List<AnnotatedRegion> matchingIntervals = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			if(value.equals(r.getAttribute(attribute))){
				matchingIntervals.add(r);
			}
		}

		return matchingIntervals;
	}

	public static List<AnnotatedRegion> SelectIntervals(GenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, String attribute, Object value){
		List<AnnotatedRegion> matchingIntervals = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end)){
			if(value.equals(r.getAttribute(attribute))){
				matchingIntervals.add(r);
			}
		}

		return matchingIntervals;
	}

	public static List<AnnotatedRegion> OverlappingIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			selected.add(r);
		}

		return selected;
	}

	public static List<AnnotatedRegion> OverlappingIntervals(GenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end)){
			selected.add(r);
		}

		return selected;
	}

	public static List<AnnotatedRegion> OverlappingIntervals(GenomicIntervalSet intervals, String chr, int start, int end){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end)){
			selected.add(r);
		}

		return selected;
	}

	public static List<AnnotatedRegion> OverlappingIntervals(StrandedGenomicIntervalSet intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			selected.add(r);
		}

		return selected;
	}

	/**
	 * whether the given interval is contained in any interval the tree
	 * 
	 * @param intervals
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public static boolean isContained(GenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end){
		boolean isContained = false;

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end)){
			isContained |= isContained(start, end, r.start, r.end);
			if(isContained)
				break;
		}

		return isContained;
	}



	/**
	 * whether the given interval is contained in any interval the tree
	 * 
	 * @param intervals
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public static boolean isContained(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		boolean isContained = false;

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			isContained |= isContained(start, end, r.start, r.end);
			if(isContained)
				break;
		}

		return isContained;
	}

	/**
	 * whether the given interval is contained in any interval the tree
	 * 
	 * @param intervals
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public static boolean isContained(GenomicIntervalSet intervals, String chr, int start, int end){
		boolean isContained = false;

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end)){
			isContained |= isContained(start, end, r.start, r.end);
			if(isContained)
				break;
		}

		return isContained;
	}

	/**
	 * whether the given interval is contained in any interval the tree
	 * 
	 * @param intervals
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	public static boolean isContained(StrandedGenomicIntervalSet intervals, String chr, int start, int end, char strand){
		boolean isContained = false;

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			isContained |= isContained(start, end, r.start, r.end);
			if(isContained)
				break;
		}

		return isContained;
	}

	public static List<AnnotatedRegion> ContainedIntervals(GenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end)){
			if(r.start >= start && r.end <= end)
				selected.add(r);
		}

		return selected;
	}

	public static List<AnnotatedRegion> ContainedIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			if(r.start >= start && r.end <= end)
				selected.add(r);
		}

		return selected;
	}
	
	public static List<AnnotatedRegion> ContainingIntervals(StrandedGenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			if(isContained(start, end, r.start, r.end))
				selected.add(r);
		}

		return selected;
	}
	
	public static List<AnnotatedRegion> ContainingIntervals(GenomicIntervalTree<Map<String,Object>> intervals, String chr, int start, int end, char strand){
		List<AnnotatedRegion> selected = new LinkedList<AnnotatedRegion>();

		for(AnnotatedRegion r : intervals.overlappingRegions(chr, start, end, strand)){
			if(isContained(start, end, r.start, r.end))
				selected.add(r);
		}

		return selected;
	}

	public static StrandedGenomicIntervalSet buildStrandedIntervalSet(Iterable<AnnotatedRegion> regions) {
		StrandedGenomicIntervalSet tree = new StrandedGenomicIntervalSet();

		for(AnnotatedRegion r : regions){			
			tree.add(r.chr, r.start, r.end, r.strand);	
		}

		return tree;
	}

	public static GenomicIntervalSet buildIntervalSet(Iterable<AnnotatedRegion> regions) {
		GenomicIntervalSet tree = new GenomicIntervalSet();

		for(AnnotatedRegion r : regions){			
			tree.add(r.chr, r.start, r.end);	
		}

		return tree;
	}

	public static void addRegions(GenomicIntervalSet tree, Iterable<AnnotatedRegion> regions) {
		for(AnnotatedRegion r : regions){			
			tree.add(r.chr, r.start, r.end);	
		}
	}

	public static void addRegions(StrandedGenomicIntervalSet tree, Iterable<AnnotatedRegion> regions) {
		for(AnnotatedRegion r : regions){			
			tree.add(r.chr, r.start, r.end, r.strand);	
		}
	}

	public static void removeRegions(GenomicIntervalSet tree, Iterable<AnnotatedRegion> regions) {
		for(AnnotatedRegion r : regions){			
			tree.remove(r.chr, r.start, r.end);	
		}
	}

	public static void removeRegions(StrandedGenomicIntervalSet tree, Iterable<AnnotatedRegion> regions) {
		for(AnnotatedRegion r : regions){			
			tree.remove(r.chr, r.start, r.end, r.strand);	
		}
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> addRegions(StrandedGenomicIntervalTree<Map<String,Object>> tree, Iterable<AnnotatedRegion> regions, boolean unique, boolean attributes) {

		for(AnnotatedRegion r : regions){
			if(!unique){
				if(attributes)
					tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
				else
					tree.add(r.chr, r.start, r.end, r.strand);
			}
			else if(!tree.contains(r.chr, r.start, r.end, r.strand)){
				if(attributes)
					tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
				else
					tree.add(r.chr, r.start, r.end, r.strand);
			}
		}

		return tree;
	}

	public static void addRegion(StrandedGenomicIntervalTree<Map<String,Object>> tree, AnnotatedRegion r, boolean unique, boolean attributes) {
		if(!unique){
			if(attributes)
				tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
			else
				tree.add(r.chr, r.start, r.end, r.strand);
		}
		else if(!tree.contains(r.chr, r.start, r.end, r.strand)){
			if(attributes)
				tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
			else
				tree.add(r.chr, r.start, r.end, r.strand);
		}
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> buildRegionsTree(Iterable<AnnotatedRegion> regions, boolean unique, boolean attributes) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			if(!unique){
				if(attributes)
					tree.add(r.chr, r.start, r.end, r.strand, r.attributes==null?new HashMap<String, Object>():r.attributes);
				else
					tree.add(r.chr, r.start, r.end, r.strand);
			}
			else if(!tree.contains(r.chr, r.start, r.end, r.strand)){
				if(attributes)
					tree.add(r.chr, r.start, r.end, r.strand, r.attributes==null?new HashMap<String, Object>():r.attributes);
				else
					tree.add(r.chr, r.start, r.end, r.strand);
			}
		}

		return tree;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> buildRegionsTree(Iterable<AnnotatedRegion> regions, boolean unique, boolean attributes, boolean annotation) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			if(!unique){
				if(attributes){
					if(annotation)
						r.addAttribute("annotation", r.annotation);
					tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
				}
				else
					tree.add(r.chr, r.start, r.end, r.strand);
			}
			else if(!tree.contains(r.chr, r.start, r.end, r.strand)){
				if(attributes){
					if(annotation)
						r.addAttribute("annotation", r.annotation);
					tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
				}
				else
					tree.add(r.chr, r.start, r.end, r.strand);
			}
		}

		return tree;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> buildRegionsTree(Iterable<AnnotatedRegion> regions, String annotation, boolean unique, boolean attributes) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			if(r.annotation.equals(annotation)){
				if(!unique){
					if(attributes)
						tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
					else
						tree.add(r.chr, r.start, r.end, r.strand);
				}
				else if(!tree.contains(r.chr, r.start, r.end, r.strand)){
					if(attributes)
						tree.add(r.chr, r.start, r.end, r.strand, r.attributes);
					else
						tree.add(r.chr, r.start, r.end, r.strand);
				}
			}
		}

		return tree;
	}

	public static GenomicIntervalTree<Map<String, Object>> buildUnstrandedRegionsTree(Iterable<AnnotatedRegion> regions, boolean unique, boolean attributes) {
		GenomicIntervalTree<Map<String,Object>> tree = new GenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			if(!unique){
				if(attributes)
					tree.add(r.chr, r.start, r.end, r.attributes);
				else
					tree.add(r.chr, r.start, r.end);
			}
			else if(!tree.contains(r.chr, r.start, r.end)){
				if(attributes)
					tree.add(r.chr, r.start, r.end, r.attributes);
				else
					tree.add(r.chr, r.start, r.end);
			}
		}

		return tree;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> buildTerminiTree(Iterable<AnnotatedRegion> regions, boolean unique, boolean attributes) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			if(!unique){
				if(attributes){
					tree.add(r.chr, r.start, r.start, r.strand, r.attributes);
					tree.add(r.chr, r.end, r.end, r.strand, r.attributes);
				}
				else{
					tree.add(r.chr, r.start, r.start, r.strand);
					tree.add(r.chr, r.end, r.end, r.strand);
				}
			}
			else if(!tree.contains(r.chr, r.start, r.end, r.strand)){
				if(attributes){
					tree.add(r.chr, r.start, r.start, r.strand, r.attributes);
					tree.add(r.chr, r.end, r.end, r.strand, r.attributes);
				}
				else{
					tree.add(r.chr, r.start, r.start, r.strand);
					tree.add(r.chr, r.end, r.end, r.strand);
				}
			}
		}

		return tree;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> buildTerminiTree(Iterable<AnnotatedRegion> regions, boolean end5p, boolean unique, boolean attributes) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			int terminus = end5p? r.get5Prime() : r.get3Prime(); 
			if(!unique){
				if(attributes){
					tree.add(r.chr, terminus, terminus, r.strand, r.attributes);
				}
				else{
					tree.add(r.chr, terminus, terminus, r.strand);
				}
			}
			else if(!tree.contains(r.chr, terminus, terminus, r.strand)){
				if(attributes){
					tree.add(r.chr, terminus, terminus, r.strand, r.attributes);
				}
				else{
					tree.add(r.chr, terminus, terminus, r.strand);
				}
			}
		}

		return tree;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> buildAttributedTerminiTree(Iterable<AnnotatedRegion> regions, boolean end5p, boolean unique) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : regions){
			int terminus = end5p? r.get5Prime() : r.get3Prime(); 
			if(!unique){
				tree.add(r.chr, terminus, terminus, r.strand, r.attributes==null?new HashMap<String, Object>():r.attributes);
			}
			else if(!tree.contains(r.chr, terminus, terminus, r.strand)){
				tree.add(r.chr, terminus, terminus, r.strand, r.attributes==null?new HashMap<String, Object>():r.attributes);
			}
		}

		return tree;
	}


	public static Comparator<AnnotatedRegion> AnnotatedRegionComparator(final boolean region5p, final boolean ascending) {

		return new Comparator<AnnotatedRegion>() {
			public int compare(AnnotatedRegion r1, AnnotatedRegion r2) {
				if(region5p){
					if(ascending)
						return r1.get5Prime() - r2.get5Prime();
					else
						return r2.get5Prime() - r1.get5Prime();
				}
				else{
					if(ascending)
						return r1.get3Prime() - r2.get3Prime();
					else
						return r2.get3Prime() - r1.get3Prime();
				}
			}
		};

	}

	public static void main(String[] args) {
		System.out.println(offsetAmount(0, 1, false, true));
		System.out.println(Util.list(offsetAmounts(1, 4, 0, 3)));
	}

	public static char opposite(char strand) {
		if(strand=='+'||strand=='-'){
			return strand=='+' ? '-' : '+';
		}
		else throw new IllegalArgumentException(Util.sprintf("No defined opposite strand for character '%c'", strand));
	}

	public static AnnotatedRegion makeRegion(String chr, int end5p, int end3p, boolean isNegativeStrand) {
		if(isNegativeStrand)
			return new AnnotatedRegion("", chr, end3p, end5p, isNegativeStrand ? '-' : '+');
		else
			return new AnnotatedRegion("", chr, end5p, end3p, isNegativeStrand ? '-' : '+');
	}

	public static int distance(int end5p, int end3p, boolean isNegativeStrand) {
		if(isNegativeStrand)
			return end5p-end3p;//new AnnotatedRegion("", chr, end3p, end5p, isNegativeStrand ? '-' : '+');
		else
			return end3p-end5p;// new AnnotatedRegion("", chr, end5p, end3p, isNegativeStrand ? '-' : '+');
	}

	public static boolean equals(AnnotatedRegion r1, AnnotatedRegion r2) {
		return r1.chr.equals(r2.chr) && r1.start == r2.start && r1.end == r2.end && r1.strand == r2.strand;
	}

	public static int terminus(int start, int end, char strand, boolean end5p) {
		Integer term = null;

		if(strand=='+'){
			if(end5p)
				term = start;
			else
				term = end;
		}
		if(strand=='-'){
			if(end5p)
				term = end;
			else
				term = start;
		}

		return term;
	}

}
