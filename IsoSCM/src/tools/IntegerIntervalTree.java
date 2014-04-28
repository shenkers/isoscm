package tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import util.Util;
/**
 * An interval tree implementation for integer data.
 * There is also a generic interval tree, but the implementation
 * was cumbersome, and included excessive object creation, so
 * hopefully this will be more efficient.
 * 
 * Implemented by extended Apache's AVL tree, tweaked to track
 * the maximum end of each interval in the subtree of each node
 * 
 * Supports the insert, delete, and overlap operations
 * 
 * It would be simple to re-implement the IntervalSet by extending
 * this class. When inserting an interval, simply find the extrema of
 * all overlaps, delete them, and insert a new overlap from one
 * extreme to the other
 * 
 * @author sol
 *
 */
public class IntegerIntervalTree<S> implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	/** Top level node. */
	private IntervalNode<S> top;

	/** Build an empty tree.
	 */
	public IntegerIntervalTree() {
		top = null;
	}

	/** Insert an element in the tree.
	 * @param element element to insert (silently ignored if null)
	 */
	public void insert(int start, int end) {
		if (top == null) {
			top = new IntervalNode<S>(start, end, null);
		} else {
			top.insert(start, end);
		}
	}
	
	/** Insert an element in the tree.
	 * @param element element to insert (silently ignored if null)
	 */
	public void insert(int start, int end, S data) {
		if (top == null) {
			top = new IntervalNode<S>(start, end, null, data);
		} else {
			top.insert(start, end, data);
		}
	}

	public List<IntervalNode<S>> overlaps(int start, int end){
		IntervalNode<S> n = top;
		if(top==null)
			return new ArrayList<IntervalNode<S>>();
		else
			return n.getOverlappingIntervals(start, end);
	}

	/** Delete an element from the tree.
	 * <p>The element is deleted only if there is a node {@code n}
	 * containing exactly the element instance specified, i.e. for which
	 * {@code n.getElement() == element}. This is purposely
	 * <em>different</em> from the specification of the
	 * {@code java.util.Set} {@code remove} method (in fact,
	 * this is the reason why a specific class has been developed).</p>
	 * @param element element to delete (silently ignored if null)
	 * @return true if the element was deleted from the tree
	 */
	public boolean delete(int start, int end) {
		for (IntervalNode<?> node = getNotSmaller(start); node != null; node = node.getNext()) {
			// loop over all elements neither smaller nor larger
			// than the specified one
			if (node.start == start && node.end == end) {
				if ((node.parent == null) && (node.left == null) && (node.right == null)) {
					// this was the last node, the tree is now empty
					//				element = null;
					top     = null;
				}
				else{
					node.delete();
				}
				return true;
			} else if (node.start > start) {
				// all the remaining elements are known to be larger,
				// the element is not in the tree
				return false;
			}
		}

		return false;
	}
	
	/** whether or not a node with the given boundaries is contained in the tree.
	 * <p>The element is deleted only if there is a node {@code n}
	 * containing exactly the element instance specified, i.e. for which
	 * {@code n.getElement() == element}. This is purposely
	 * <em>different</em> from the specification of the
	 * {@code java.util.Set} {@code remove} method (in fact,
	 * this is the reason why a specific class has been developed).</p>
	 * @param element element to delete (silently ignored if null)
	 * @return true if the element was deleted from the tree
	 */
	public boolean contains(int start, int end) {
		for (IntervalNode<?> node = getNotSmaller(start); node != null; node = node.getNext()) {
			// loop over all elements neither smaller nor larger
			// than the specified one
			if (node.start == start && node.end == end) {
				return true;
			} else if (node.start > start) {
				// all the remaining elements are known to be larger,
				// the element is not in the tree
				return false;
			}
		}

		return false;
	}
	
	/** Returns a node with the given boundaries if it is contained in the tree, and null otherwise.
	 * <p>The element is deleted only if there is a node {@code n}
	 * containing exactly the element instance specified, i.e. for which
	 * {@code n.getElement() == element}. This is purposely
	 * <em>different</em> from the specification of the
	 * {@code java.util.Set} {@code remove} method (in fact,
	 * this is the reason why a specific class has been developed).</p>
	 * @param element element to delete (silently ignored if null)
	 * @return true if the element was deleted from the tree
	 */
	public IntervalNode<S> get(int start, int end) {
		for (IntervalNode<S> node = getNotSmaller(start); node != null; node = node.getNext()) {
			// loop over all elements neither smaller nor larger
			// than the specified one
			if (node.start == start && node.end == end) {
				return node;
			} else if (node.start > start) {
				// all the remaining elements are known to be larger,
				// the element is not in the tree
				return null;
			}
		}

		return null;
	}

	/** Check if the tree is empty.
	 * @return true if the tree is empty
	 */
	public boolean isEmpty() {
		return top == null;
	}


	/** Get the number of elements of the tree.
	 * @return number of elements contained in the tree
	 */
	public int size() {
		return (top == null) ? 0 : top.size();
	}

	/** Get the node whose element is the smallest one in the tree.
	 * @return the tree node containing the smallest element in the tree
	 * or null if the tree is empty
	 * @see #getLargest
	 * @see #getNotSmaller
	 * @see #getNotLarger
	 * @see Node#getPrevious
	 * @see Node#getNext
	 */
	public IntervalNode<S> getSmallest() {
		return (top == null) ? null : top.getSmallest();
	}

	/** Get the node whose element is the largest one in the tree.
	 * @return the tree node containing the largest element in the tree
	 * or null if the tree is empty
	 * @see #getSmallest
	 * @see #getNotSmaller
	 * @see #getNotLarger
	 * @see Node#getPrevious
	 * @see Node#getNext
	 */
	public IntervalNode<S> getLargest() {
		return (top == null) ? null : top.getLargest();
	}

	/** Get the node whose element is not smaller than the reference object.
	 * @param reference reference object (may not be in the tree)
	 * @return the tree node containing the smallest element not smaller
	 * than the reference object or null if either the tree is empty or
	 * all its elements are smaller than the reference object
	 * @see #getSmallest
	 * @see #getLargest
	 * @see #getNotLarger
	 * @see Node#getPrevious
	 * @see Node#getNext
	 */
	public IntervalNode<S> getNotSmaller(int start) {
		IntervalNode<S> candidate = null;
		for (IntervalNode<S> node = top; node != null;) {
			if (node.start < start) {
				if (node.right == null) {
					return candidate;
				}
				node = node.right;
			} else {
				candidate = node;
				if (node.left == null) {
					return candidate;
				}
				node = node.left;
			}
		}
		return null;
	}

	/** Get the node whose element is not larger than the reference object.
	 * @param reference reference object (may not be in the tree)
	 * @return the tree node containing the largest element not larger
	 * than the reference object (in which case the node is guaranteed
	 * not to be empty) or null if either the tree is empty or all its
	 * elements are larger than the reference object
	 * @see #getSmallest
	 * @see #getLargest
	 * @see #getNotSmaller
	 * @see Node#getPrevious
	 * @see Node#getNext
	 */
	public IntervalNode<S> getNotLarger(int start) {
		IntervalNode<S> candidate = null;
		for (IntervalNode<S> node = top; node != null;) {
			if (node.start > start) {
				if (node.left == null) {
					return candidate;
				}
				node = node.left;
			} else {
				candidate = node;
				if (node.right == null) {
					return candidate;
				}
				node = node.right;
			}
		}
		return null;
	}

	public IntervalNode<S> getNotLargerEnd(int end) {
		IntervalNode<S> candidate = null;
		for (IntervalNode<S> node = top; node != null;) {
			if (node.start > end) {
				if (node.left == null) {
					return candidate;
				}
				node = node.left;
			} else {
				candidate = node;
				if (node.right == null) {
					return candidate;
				}
				node = node.right;
			}
		}
		return null;
	}
	
	
	public String toNewickString() {
		return top.toNewickString();
	}
	
	public static void main(String[] args) {
		IntegerIntervalTree<Integer> avlt = new IntegerIntervalTree<Integer>();
		for(int i=0; i<7; i++){
			int start = (int) (Math.random()*20);
			int end = (int) (Math.random()*5) + 1;
			avlt.insert(start, start+end);
			System.out.println("inserting:"+Util.list(start,start+end));
		}
		System.out.println();
		for(int i=0; i<5; i++){
			int start = (int) (Math.random()*20);
			int end = (int) (Math.random()*5) + 1;
			System.out.println("overlapping:"+Util.list(start,start+end));
			System.out.println(avlt.overlaps(start, start+end));
//			avlt.overlaps(start, start+end);
					System.out.println();
		}
		avlt.overlaps(2, 2);
		
		System.out.println(avlt.toNewickString());
		for(IntervalNode<Integer> n : avlt.overlaps(2, 2)){
			System.out.println("deleting:"+n);
			avlt.delete(n.start, n.end);
		}
		System.out.println(avlt.toNewickString());
	}

}


