package tools;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import util.Util;


public class IntegerIntervalSet<S> implements Serializable, Iterable<int[]>{

	private static final long serialVersionUID = 1L;

	IntegerIntervalTree<S> tree;

	public IntegerIntervalSet() {
		tree = new IntegerIntervalTree<S>();
	}

	public void add(int start, int end){
		if (end<start) {
			throw new IllegalArgumentException(Util.sprintf("[%d,%d] is an invalid interval", start,end));
		}

		int minStart = start; int maxEnd = end;

		List<IntervalNode<S>> overlaps = tree.overlaps(start-1, end+1);

		List<Integer[]> toDelete = new LinkedList<Integer[]>();
		for(IntervalNode<S> overlap : overlaps){
			minStart = Math.min(overlap.start,minStart);
			maxEnd = Math.max(overlap.end,maxEnd);
			toDelete.add(new Integer[]{overlap.start, overlap.end});
		}
		
		for(Integer[] overlapped : toDelete){
			tree.delete(overlapped[0],overlapped[1]);
		}

		tree.insert(minStart, maxEnd);	
	}

	// TODO could add facility for generic associated data, but would need a merge method
	// when an insert combines two regions 

	public void remove(int start, int end){
		if (end<start) {
			throw new IllegalArgumentException();
		}

		List<IntervalNode<S>> overlaps = tree.overlaps(start, end);
		List<int[]> toAdd = new LinkedList<int[]>();
		
		int[][] to_delete = new int[overlaps.size()][2];

		for(int i=0; i<overlaps.size(); i++){
			IntervalNode<S> overlap = overlaps.get(i); 
			int overlap_start = overlap.start;
			int overlap_end = overlap.end;
			to_delete[i][0] = overlap_start;
			to_delete[i][1] = overlap_end;

			if(overlap_start<start){
//				tree.insert(overlap_start, start-1);
				toAdd.add(new int[]{overlap_start, start-1});
			}
			if(overlap_end>end){
//				tree.insert(end+1, overlap_end);
				toAdd.add(new int[]{end+1, overlap_end});
			}
		}

		for(int i=0; i<to_delete.length; i++){
			tree.delete(to_delete[i][0], to_delete[i][1]);		
		}
		for(int[] interval : toAdd){
			tree.insert(interval[0], interval[1]);
		}
	}

	public Iterator<int[]> iterator() {
		return new IntervalIterator(tree.getSmallest());
	}

	class IntervalIterator implements Iterator<int[]>{
		IntervalNode<S> node;

		public IntervalIterator(IntervalNode<S> smallest) {
			node = smallest;
		}

		public boolean hasNext() {
			return node!=null;
		}

		public int[] next() {
			int[] interval = new int[]{node.start,node.end};
			node = node.getNext();
			return interval;
		}

		public void remove() {
			// TODO Auto-generated method stub

		}

	}

	public String toString() {
		List<List<Integer>> intervals = new ArrayList<List<Integer>>(tree.size());

		IntervalNode<S> node = tree.getSmallest();
		while(node!=null){	
			intervals.add(Util.list(node.start,node.end));
			node = node.getNext();
		}

		return intervals.toString();
	}

	public int size(){
		return tree.size();
	}

	/**
	 * 
	 * @param start
	 * @param end
	 * @return the set of intervals that overlap with the given segment
	 */
	public List<List<Integer>> overlappingIntervals(int start, int end){
		List<List<Integer>> overlappingIntervals = new LinkedList<List<Integer>>();

		for(IntervalNode<S> overlap : tree.overlaps(start, end)){
			overlappingIntervals.add(Util.list(overlap.start,overlap.end));	
		}

		return overlappingIntervals;
	}

	public Iterable<List<Integer>> overlappingIntervalIterable(int start, int end){
		return overlappingIntervals(start, end);
	}

	public static void main(String[] args) {
		List<Integer> li = Util.list(1,3,4);
		System.out.println(Collections.binarySearch(li, 0));

		IntegerIntervalSet is = new IntegerIntervalSet();
		IntegerIntervalSet is2 = new IntegerIntervalSet();

		is2.add(10,13);
		is2.add(16, 20);

		is.add(-4, -3);
		System.out.println(is);
		is.add(0, 11);
		System.out.println(is);
		is.add(-1, 10);
		System.out.println(is);
		is.add(-2, -2);
		System.out.println(is);
		is.remove(-2, -2);
		System.out.println(is);
		is.remove(-5, -4);
		System.out.println(is);

		//		System.out.println(is);
		is.add(15, 16);
		System.out.println(is);
		System.out.println(is2);
		System.out.println("!");

	}


}
