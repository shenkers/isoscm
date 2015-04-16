package tools;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import util.Util;
import util.Util.ExtremeTracker;

public class SegmentFragmenter {

	/**
	 * assumes that segments are all added to interval tree through this method.
	 * otherwise it may reach an inconsistent state and behavior is unpredictable.
	 * 
	 * @param itervalTree
	 * @param chr
	 * @param start
	 * @param end
	 */
	public static void add(GenomicIntervalTree<Map<String,Object>> itervalTree, String chr, int start, int end){

		List<AnnotatedRegion> overlaps = new LinkedList<AnnotatedRegion>();
	Set<Integer> startBreakPoints = new HashSet<Integer>();
		Set<Integer> endBreakPoints = new HashSet<Integer>();
		startBreakPoints.add(start);
		endBreakPoints.add(start-1);
		endBreakPoints.add(end);
		startBreakPoints.add(end+1);
		ExtremeTracker<Integer> extremes = new ExtremeTracker<Integer>(new Util.ComparableComparator<Integer>());
		for(AnnotatedRegion r : itervalTree.overlappingRegions(chr, start, end)){
			overlaps.add(r);
			extremes.put(r.start);
			extremes.put(r.end);
			endBreakPoints.add(r.start-1);
			startBreakPoints.add(r.start);
			endBreakPoints.add(r.end);
			startBreakPoints.add(r.end+1);
		}

		for(AnnotatedRegion r : overlaps){
			itervalTree.remove(r.chr, r.start, r.end);
		}

		List<Integer> startBreaks = Util.list(startBreakPoints);
		List<Integer> endBreaks = Util.list(endBreakPoints);
		Collections.sort(startBreaks);
		Collections.sort(endBreaks);

		for(int i=0; i<startBreaks.size()-1; i+=1){
//			System.out.println(Util.list(startBreaks.get(i),endBreaks.get(i+1)));
			itervalTree.add(chr, startBreaks.get(i),endBreaks.get(i+1));
		}	
	}

	public static void add(StrandedGenomicIntervalTree<Map<String,Object>> intervalTree, String chr, int start, int end, char strand){

		List<AnnotatedRegion> overlaps = new LinkedList<AnnotatedRegion>();
		Set<Integer> startBreakPoints = new HashSet<Integer>();
		Set<Integer> endBreakPoints = new HashSet<Integer>();
		startBreakPoints.add(start);
		endBreakPoints.add(start-1);
		endBreakPoints.add(end);
		startBreakPoints.add(end+1);
		ExtremeTracker<Integer> extremes = new ExtremeTracker<Integer>(new Util.ComparableComparator<Integer>());
		for(AnnotatedRegion r : intervalTree.overlappingRegions(chr, start, end, strand)){
			overlaps.add(r);
			extremes.put(r.start);
			extremes.put(r.end);
			endBreakPoints.add(r.start-1);
			startBreakPoints.add(r.start);
			endBreakPoints.add(r.end);
			startBreakPoints.add(r.end+1);
		}

		for(AnnotatedRegion r : overlaps){
			intervalTree.remove(r.chr, r.start, r.end, strand);
		}

		List<Integer> startBreaks = Util.list(startBreakPoints);
		List<Integer> endBreaks = Util.list(endBreakPoints);
		Collections.sort(startBreaks);
		Collections.sort(endBreaks);



		for(int i=0; i<startBreaks.size()-1; i+=1){
//			System.out.println(Util.list(startBreaks.get(i),endBreaks.get(i+1)));
			intervalTree.add(chr, startBreaks.get(i),endBreaks.get(i+1), strand);
		}		
	}
	
	public static void fragment(GenomicIntervalTree<Map<String,Object>> intervalTree, String chr, int pos, boolean break_before){

		List<AnnotatedRegion> overlaps = new LinkedList<AnnotatedRegion>();
		Set<Integer> startBreakPoints = new HashSet<Integer>();
		Set<Integer> endBreakPoints = new HashSet<Integer>();
	
		if(break_before){
			startBreakPoints.add(pos);
			endBreakPoints.add(pos-1);
			}
		else{
			startBreakPoints.add(pos+1);
			endBreakPoints.add(pos);
			}
		ExtremeTracker<Integer> extremes = new ExtremeTracker<Integer>(new Util.ComparableComparator<Integer>());
		for(AnnotatedRegion r : intervalTree.overlappingRegions(chr, pos, pos)){
			overlaps.add(r);
			extremes.put(r.start);
			extremes.put(r.end);
			endBreakPoints.add(r.start-1);
			startBreakPoints.add(r.start);
			endBreakPoints.add(r.end);
			startBreakPoints.add(r.end+1);
		}

		for(AnnotatedRegion r : overlaps){
			intervalTree.remove(r.chr, r.start, r.end);
		}

		List<Integer> startBreaks = Util.list(startBreakPoints);
		List<Integer> endBreaks = Util.list(endBreakPoints);
		Collections.sort(startBreaks);
		Collections.sort(endBreaks);



		for(int i=0; i<startBreaks.size()-1; i+=1){
//			System.out.println(Util.list(startBreaks.get(i),endBreaks.get(i+1)));
			intervalTree.add(chr, startBreaks.get(i),endBreaks.get(i+1));
		}		
	}
	
	public static void fragment(StrandedGenomicIntervalTree<Map<String,Object>> intervalTree, String chr, int pos, char strand, boolean break_before){

		List<AnnotatedRegion> overlaps = new LinkedList<AnnotatedRegion>();
		Set<Integer> startBreakPoints = new HashSet<Integer>();
		Set<Integer> endBreakPoints = new HashSet<Integer>();
	
		if(break_before){
			startBreakPoints.add(pos);
			endBreakPoints.add(pos-1);
			}
		else{
			startBreakPoints.add(pos+1);
			endBreakPoints.add(pos);
			}
		ExtremeTracker<Integer> extremes = new ExtremeTracker<Integer>(new Util.ComparableComparator<Integer>());
		for(AnnotatedRegion r : intervalTree.overlappingRegions(chr, pos, pos, strand)){
			overlaps.add(r);
			extremes.put(r.start);
			extremes.put(r.end);
			endBreakPoints.add(r.start-1);
			startBreakPoints.add(r.start);
			endBreakPoints.add(r.end);
			startBreakPoints.add(r.end+1);
		}

		for(AnnotatedRegion r : overlaps){
			intervalTree.remove(r.chr, r.start, r.end, strand);
		}

		List<Integer> startBreaks = Util.list(startBreakPoints);
		List<Integer> endBreaks = Util.list(endBreakPoints);
		Collections.sort(startBreaks);
		Collections.sort(endBreaks);



		for(int i=0; i<startBreaks.size()-1; i+=1){
//			System.out.println(Util.list(startBreaks.get(i),endBreaks.get(i+1)));
			intervalTree.add(chr, startBreaks.get(i),endBreaks.get(i+1), strand);
		}		
	}

	public static void main(String[] args) {
		StrandedGenomicIntervalTree<Map<String,Object>> tree = new StrandedGenomicIntervalTree<Map<String,Object>>();

		tree.add("1", 1, 10,'.');
		//		
				tree.add("1", 15, 20,'.');
		//		
		//		add(tree, "1", 8, 17);
		add(tree, "1", 10, 17,'.');
		add(tree, "1", 10, 17,'.');
		fragment(tree, "1", 12, '.', false);
		//		System.out.println(tree);
		for(AnnotatedRegion r : tree)
			System.out.println(r);
		
		System.out.println(tree.size());
	}

}
