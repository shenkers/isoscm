package tools;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.geometry.partitioning.utilities.AVLTree;

/** This class implements AVL trees nodes.
 * <p>AVL tree nodes implement all the logical structure of the
 * tree. Nodes are created by the {@link AVLTree AVLTree} class.</p>
 * <p>The nodes are not independant from each other but must obey
 * specific balancing constraints and the tree structure is
 * rearranged as elements are inserted or deleted from the tree. The
 * creation, modification and tree-related navigation methods have
 * therefore restricted access. Only the order-related navigation,
 * reading and delete methods are public.</p>
 * @see AVLTree
 */
public class IntervalNode<S> implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/** Enum for tree skew factor. */
	private static enum Skew {
		/** Code for left high trees. */
		LEFT_HIGH,

		/** Code for right high trees. */
		RIGHT_HIGH,

		/** Code for Skew.BALANCED trees. */
		BALANCED;
	}

	/** Left sub-tree. */
	IntervalNode<S> left;

	/** Right sub-tree. */
	IntervalNode<S> right;

	/** Parent tree. */
	IntervalNode<S> parent;

	/** Skew factor. */
	Skew skew;

	int maxEnd;

	int start;

	int end;

	boolean rebalance;

	public S data;

	/** Build a node for a specified element.
	 * @param element element
	 * @param parent parent node
	 */
	IntervalNode(int start, int end, final IntervalNode<S> parent) {
		left         = null;
		right        = null;
		rebalance=false;
		this.parent  = parent;
		skew         = Skew.BALANCED;
		maxEnd = end;
		this.start=start;
		this.end=end;
	}
	
	/** Build a node for a specified element.
	 * @param element element
	 * @param parent parent node
	 */
	IntervalNode(int start, int end, final IntervalNode<S> parent, S data) {
		left         = null;
		right        = null;
		rebalance=false;
		this.parent  = parent;
		skew         = Skew.BALANCED;
		maxEnd = end;
		this.start=start;
		this.end=end;
		this.data=data;
	}

	public int getStart(){
		return start;
	}

	public int getEnd(){
		return end;
	}

	boolean overlaps(int start, int end){
		return (this.start >= start && this.start <= end) ||
				(this.end >= start && this.end <= end) ||
				(this.start <= start && this.end >= end);
	}

	public List<IntervalNode<S>> getOverlappingIntervals(int start, int end) {
		List<IntervalNode<S>> overlaps = new LinkedList<IntervalNode<S>>();
		return getOverlappingIntervals(start, end, overlaps);
	}

	public List<IntervalNode<S>> getOverlappingIntervals(int start, int end, List<IntervalNode<S>> overlaps) {
		if(overlaps(start,end)){
			overlaps.add(this);
		}

		if(left != null && left.maxEnd >= start){
			left.getOverlappingIntervals(start,end,overlaps);
		}
		if(this.start <= end && right != null && right.maxEnd >= start){
			right.getOverlappingIntervals(start,end,overlaps);
		}

		return overlaps;
	}

	/** Get the number of elements of the tree rooted at this node.
	 * @return number of elements contained in the tree rooted at this node
	 */
	int size() {
		return 1 + ((left  == null) ? 0 : left.size()) + ((right == null) ? 0 : right.size());
	}

	/** Get the node whose element is the smallest one in the tree
	 * rooted at this node.
	 * @return the tree node containing the smallest element in the
	 * tree rooted at this node or null if the tree is empty
	 * @see #getLargest
	 */
	IntervalNode<S> getSmallest() {
		IntervalNode<S> node = this;
		while (node.left != null) {
			node = node.left;
		}
		return node;
	}

	/** Get the node whose element is the largest one in the tree
	 * rooted at this node.
	 * @return the tree node containing the largest element in the
	 * tree rooted at this node or null if the tree is empty
	 * @see #getSmallest
	 */
	IntervalNode<S> getLargest() {
		IntervalNode<S> node = this;
		while (node.right != null) {
			node = node.right;
		}
		return node;
	}

	/** Get the node containing the next smaller or equal element.
	 * @return node containing the next smaller or equal element or
	 * null if there is no smaller or equal element in the tree
	 * @see #getNext
	 */
	public IntervalNode<S> getPrevious() {

		if (left != null) {
			final IntervalNode<S> node = left.getLargest();
			if (node != null) {
				return node;
			}
		}

		for (IntervalNode<S> node = this; node.parent != null; node = node.parent) {
			if (node != node.parent.left) {
				return node.parent;
			}
		}

		return null;

	}

	/** Get the node containing the next larger or equal element.
	 * @return node containing the next larger or equal element (in
	 * which case the node is guaranteed not to be empty) or null if
	 * there is no larger or equal element in the tree
	 * @see #getPrevious
	 */
	public IntervalNode<S> getNext() {

		if (right != null) {
			final IntervalNode<S> node = right.getSmallest();
			if (node != null) {
				return node;
			}
		}

		for (IntervalNode<S> node = this; node.parent != null; node = node.parent) {
			if (node != node.parent.right) {
				return node.parent;
			}
		}

		return null;

	}

	/** Insert an element in a sub-tree.
	 * @param newElement element to insert
	 * @return true if the parent tree should be re-Skew.BALANCED
	 */
	IntervalNode<S> insert(int start, int end) {
		if(this.maxEnd < end)
			maxEnd = end;
		
		if (start < this.start) {
			// the inserted element is smaller than the node
			if (left == null) {
				left = new IntervalNode<S>(start, end, this);
				rebalance = rebalanceLeftGrown();
				return left;
			}
			IntervalNode<S> node = left.insert(start, end);
			if(left.rebalance)
				rebalance = rebalanceLeftGrown();
			else
				rebalance = false;
			return node;
		}

		// the inserted element is equal to or greater than the node
	
		if (right == null) {
			right = new IntervalNode<S>(start,end, this);
			rebalance = rebalanceRightGrown();
			return right;
		}
		IntervalNode<S> node = right.insert(start, end);
		if(right.rebalance)
			rebalance = rebalanceRightGrown();
		else
			rebalance = false;
		return node;

	}
	
	IntervalNode<S> insert(int start, int end, S data) {
		if(this.maxEnd < end)
			maxEnd = end;
		
		if (start < this.start) {
			// the inserted element is smaller than the node
			if (left == null) {
				left = new IntervalNode<S>(start, end, this, data);
//				System.out.println("rebalanc left 1");
				rebalance = rebalanceLeftGrown();
				return left;
			}
			IntervalNode<S> node = left.insert(start, end, data);
			if(left.rebalance){
//				System.out.println("rebalanc left 2");
					rebalance = rebalanceLeftGrown();
			}else
				rebalance = false;
			return node;
		}

		// the inserted element is equal to or greater than the node
		if(this.maxEnd < end)
			maxEnd = end;
		if (right == null) {
			right = new IntervalNode<S>(start,end, this, data);
//			System.out.println("rebalanc right 1");
			rebalance = rebalanceRightGrown();
			return right;
		}
		IntervalNode<S> node = right.insert(start, end, data);
		if(right.rebalance){
//			System.out.println("rebalanc right 2");
				rebalance = rebalanceRightGrown();
		}else
			rebalance = false;
		return node;

	}

	/** Delete the node from the tree.
	 */
	public void delete() {
		IntervalNode<S> node;
		IntervalNode<S> child;
		boolean leftShrunk;
		if ((left == null) && (right == null)) {
			node       = this;
			//					element    = null;
			leftShrunk = node == node.parent.left;
			child      = null;
		} else {
			node       = (left != null) ? left.getLargest() : right.getSmallest();
			this.start = node.start;
			this.end = node.end;
			this.data = node.data;
			leftShrunk = node == node.parent.left;
			child      = (node.left != null) ? node.left : node.right;
		}

		node = node.parent;
		if (leftShrunk) {
			node.left = child;
		} else {
			node.right = child;
		}
		if (child != null) {
			child.parent = node;
		}

		node.propogateUpdate();
		this.propogateUpdate();
		
		while (leftShrunk ? node.rebalanceLeftShrunk() : node.rebalanceRightShrunk()) {
			if (node.parent == null) {
				return;
			}
			leftShrunk = node == node.parent.left;
			node = node.parent;
		}
	}

	/** Re-balance the instance as left sub-tree has grown.
	 * @return true if the parent tree should be reSkew.BALANCED too
	 */
	private boolean rebalanceLeftGrown() {
		switch (skew) {
		case LEFT_HIGH:
			if (left.skew == Skew.LEFT_HIGH) {
				rotateCW();
				skew       = Skew.BALANCED;
				right.skew = Skew.BALANCED;
			} else {
				final Skew s = left.right.skew;
				left.rotateCCW();
				rotateCW();
				switch(s) {
				case LEFT_HIGH:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.RIGHT_HIGH;
					break;
				case RIGHT_HIGH:
					left.skew  = Skew.LEFT_HIGH;
					right.skew = Skew.BALANCED;
					break;
				default:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.BALANCED;
				}
				skew = Skew.BALANCED;
			}
			return false;
		case RIGHT_HIGH:
			skew = Skew.BALANCED;
			return false;
		default:
			skew = Skew.LEFT_HIGH;
			return true;
		}
	}

	/** Re-balance the instance as right sub-tree has grown.
	 * @return true if the parent tree should be reSkew.BALANCED too
	 */
	private boolean rebalanceRightGrown() {
		switch (skew) {
		case LEFT_HIGH:
			skew = Skew.BALANCED;
			return false;
		case RIGHT_HIGH:
			if (right.skew == Skew.RIGHT_HIGH) {
				rotateCCW();
				skew      = Skew.BALANCED;
				left.skew = Skew.BALANCED;
			} else {
				final Skew s = right.left.skew;
				right.rotateCW();
				rotateCCW();
				switch (s) {
				case LEFT_HIGH:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.RIGHT_HIGH;
					break;
				case RIGHT_HIGH:
					left.skew  = Skew.LEFT_HIGH;
					right.skew = Skew.BALANCED;
					break;
				default:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.BALANCED;
				}
				skew = Skew.BALANCED;
			}
			return false;
		default:
			skew = Skew.RIGHT_HIGH;
			return true;
		}
	}

	/** Re-balance the instance as left sub-tree has shrunk.
	 * @return true if the parent tree should be reSkew.BALANCED too
	 */
	private boolean rebalanceLeftShrunk() {
		switch (skew) {
		case LEFT_HIGH:
			skew = Skew.BALANCED;
			return true;
		case RIGHT_HIGH:
			if (right.skew == Skew.RIGHT_HIGH) {
				rotateCCW();
				skew      = Skew.BALANCED;
				left.skew = Skew.BALANCED;
				return true;
			} else if (right.skew == Skew.BALANCED) {
				rotateCCW();
				skew      = Skew.LEFT_HIGH;
				left.skew = Skew.RIGHT_HIGH;
				return false;
			} else {
				final Skew s = right.left.skew;
				right.rotateCW();
				rotateCCW();
				switch (s) {
				case LEFT_HIGH:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.RIGHT_HIGH;
					break;
				case RIGHT_HIGH:
					left.skew  = Skew.LEFT_HIGH;
					right.skew = Skew.BALANCED;
					break;
				default:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.BALANCED;
				}
				skew = Skew.BALANCED;
				return true;
			}
		default:
			skew = Skew.RIGHT_HIGH;
			return false;
		}
	}

	/** Re-balance the instance as right sub-tree has shrunk.
	 * @return true if the parent tree should be reSkew.BALANCED too
	 */
	private boolean rebalanceRightShrunk() {
		switch (skew) {
		case RIGHT_HIGH:
			skew = Skew.BALANCED;
			return true;
		case LEFT_HIGH:
			if (left.skew == Skew.LEFT_HIGH) {
				rotateCW();
				skew       = Skew.BALANCED;
				right.skew = Skew.BALANCED;
				return true;
			} else if (left.skew == Skew.BALANCED) {
				rotateCW();
				skew       = Skew.RIGHT_HIGH;
				right.skew = Skew.LEFT_HIGH;
				return false;
			} else {
				final Skew s = left.right.skew;
				left.rotateCCW();
				rotateCW();
				switch (s) {
				case LEFT_HIGH:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.RIGHT_HIGH;
					break;
				case RIGHT_HIGH:
					left.skew  = Skew.LEFT_HIGH;
					right.skew = Skew.BALANCED;
					break;
				default:
					left.skew  = Skew.BALANCED;
					right.skew = Skew.BALANCED;
				}
				skew = Skew.BALANCED;
				return true;
			}
		default:
			skew = Skew.LEFT_HIGH;
			return false;
		}
	}

	/** Perform a clockwise rotation rooted at the instance.
	 * <p>The skew factor are not updated by this method, they
	 * <em>must</em> be updated by the caller</p>
	 */
	private void rotateCW() {

		// swap this with left
		int tStart = start; 
		int tEnd = end;
		S tData = data;
		
		start=left.start; 
		end = left.end;
		data = left.data;
		
		left.start = tStart; 
		left.end = tEnd;
		left.data = tData;

		// adjust left children
		final IntervalNode<S> tmpNode   = left;
		left                 = tmpNode.left;
		tmpNode.left         = tmpNode.right;
		tmpNode.right        = right;
		right                = tmpNode;

		if (left != null) {
			left.parent = this;
		}
		if (right.right != null) {
			right.right.parent = right;
		}

		if(right!=null)
			right.updateMaxEnd();
		updateMaxEnd();			

	}

	/** Perform a counter-clockwise rotation rooted at the instance.
	 * <p>The skew factor are not updated by this method, they
	 * <em>must</em> be updated by the caller</p>
	 */
	private void rotateCCW() {

		int tStart = start; 
		int tEnd = end;
		S tData = data;
		
		start=right.start; 
		end = right.end;
		data = right.data;
		
		right.start = tStart; 
		right.end = tEnd;
		right.data = tData;
		
		final IntervalNode<S> tmpNode    = right;
		right                 = tmpNode.right;
		tmpNode.right         = tmpNode.left;
		tmpNode.left          = left;
		left                  = tmpNode;

		if (right != null) {
			right.parent = this;
		}
		if (left.left != null) {
			left.left.parent = left;
		}

		if(left!=null)
			left.updateMaxEnd();
		updateMaxEnd();

	}

	private boolean updateMaxEnd(){
		int max = end;
		if(left != null && max < left.maxEnd)
			max = left.maxEnd;
		if(right != null && max < right.maxEnd)
			max = right.maxEnd;

		if(max != maxEnd){
			maxEnd = max;
			return true;
		}
		else{
			return false;
		}

	}

	private void propogateUpdate(){
		boolean updated = updateMaxEnd();
		if(updated && parent != null){
			parent.propogateUpdate();
		}
	}

	public String toString(){
		return "["+start+", "+end+"]";
	}

	public String toNewickString() {
		/*
	        We're going to assume a pre-order traversal by default
		 */

		String stringRepresentation = "("+toString()+",";

		if(left != null) {
			stringRepresentation += left.toNewickString();
		}
		stringRepresentation += ",";
		if(right != null) {
			stringRepresentation += right.toNewickString();
		}

		stringRepresentation +=  ")";

		return stringRepresentation;
	}

}