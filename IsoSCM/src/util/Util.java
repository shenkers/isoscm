package util;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.util.ArithmeticUtils;


public class Util {
	
	public static class CriteriaIterable<T> implements Iterable<T>{

		Iterable<T> iterable;
		Criteria<T> c;

		public CriteriaIterable(Iterable<T> iterable, Criteria<T> c){
			this.iterable=iterable;
			this.c=c;
		}

		public Iterator<T> iterator() {
			return new CriteriaIterator<T>(iterable.iterator(), c);
		}

	}

	public static class CriteriaIterator<T> implements Iterator<T>{

		Iterator<T> it;
		T next;
		boolean hasNext;
		Criteria<T> c;

		public CriteriaIterator(	Iterator<T> it, Criteria<T> c) {
			this.it = it;
			this.c=c;
			findNext();
		}

		public boolean hasNext() {
			return hasNext;
		}

		void findNext(){
			while(it.hasNext()){
				T n = it.next();
				if(c.satisfies(n)){
					next=n;
					hasNext=true;
					return;
				}
			}
			next=null;
			hasNext=false;
		}

		public T next() {
			T last = next;
			findNext();
			return last;
		}

		public void remove() {
			// TODO Auto-generated method stub
		}

	}

	public static class CircularArray<T>{	
		T[] data;
		int i;
		int l;
		
		public CircularArray(int l, Class<T> c){
			data = (T[]) Array.newInstance(c, l);
			i=0;
			this.l = l;
		}
		
		public void set(int index, T value){
			if(index < 0 || index > l-1)
				throw new IndexOutOfBoundsException(Util.sprintf("index '%d' is not in range [0,%d]",index,l-1));
			data[(index+i)%l] = value;
		}
		
		public T get(int index){
			if(index < 0 || index > l-1)
				throw new IndexOutOfBoundsException(Util.sprintf("index '%d' is not in range [0,%d]",index,l-1));
			return data[(index+i)%l];
		}
		
		public void shift_left(){
			data[i] = null;
			i = (i+1)%l;
		}
		
		public void shift_right(){
			i = (i+l-1)%l;
			data[i] = null;
		}
		
		public void rotate_left(){
			i = (i+1)%l;
		}
		
		public void rotate_right(){
			i = (i+l-1)%l;
		}
		
		public T[] linearize(){
			Class<? extends T[]> c = (Class<? extends T[]>) data.getClass();
			T[] d = (T[]) Array.newInstance(c.getComponentType(), l);
			for(int j=0; j<l; j++){
				d[j] = data[(i+j)%l];
			}
			return d;
		}
	}
	
	public static class IntegralTable{
		
		double[] integral;
		
		public IntegralTable(double[] y) {
			integral = new double[y.length+1];
			for(int i=1; i<y.length+1; i++){
				integral[i] = integral[i-1]+y[i-1];
			}
		}
		
		/**
		 * 
		 * @param i
		 * @param j
		 * @return the sum of y from index i to j
		 */
		public double sum(int i, int j){
			return integral[j+1] - integral[i];
		}
	}

	public static TreeNode parseNewick(String newickString, TreeNode root){
		int bracket = 0;
		int from = 0;
		int closed = -1;

		for(int leftPointer = 0; leftPointer < newickString.length(); leftPointer++){
			if(newickString.charAt(leftPointer) == '(')
				bracket++;
			else if(newickString.charAt(leftPointer) == ')'){
				bracket--;
				closed = leftPointer;
			}

			if(bracket == 0 && newickString.charAt(leftPointer) == ',' ){
				if(closed != -1)					
					root.addChild(parseNewick(newickString.substring(from + 1, closed), new TreeNode(newickString.substring(closed + 1,leftPointer))));			
				else
					root.addChild(new TreeNode(newickString.substring(from, leftPointer)));

				from = leftPointer + 1;
				closed = -1;
			}
			if(leftPointer == newickString.length()-1){
				if(closed != -1)	
					root.addChild(parseNewick(newickString.substring(from + 1, closed), new TreeNode(newickString.substring(closed + 1,leftPointer+1))));

				else
					root.addChild(new TreeNode(newickString.substring(from, leftPointer+1)));

				from = leftPointer + 1;
				closed = -1;
			}
		}


		return root;		
	}


	public static class OrderedMap<K, V> extends TreeMap<K, V> implements Serializable{
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		AdditionOrderComparator<K> additionComparator;
		public OrderedMap() {
			super(new AdditionOrderComparator<K>());
			additionComparator=(AdditionOrderComparator<K>) this.comparator();
		}

		public V put(K key, V value) {
			super.remove(key);
			//			additionComparator.remove(key);
			additionComparator.add(key);
			return super.put(key, value);
		};

	}

	private static class AdditionOrderComparator<K> implements Comparator<K>, Serializable{

		Map<K,Integer> additionOrder;
		int i;

		public AdditionOrderComparator() {
			additionOrder=new HashMap<K, Integer>();
			i=1;
		}

		public int compare(K o1, K o2) {
			if(additionOrder.containsKey(o1) && additionOrder.containsKey(o2))
				return additionOrder.get(o1)-additionOrder.get(o2);
			else
				if(!(additionOrder.containsKey(o1) || additionOrder.containsKey(o2)))
					return 0;
				else
					return additionOrder.containsKey(o1) ? 1 : -1;
		}

		public void add(K k){
			additionOrder.put(k, i);
			i++;
		}

		//		public void remove(K k){
		//			additionOrder.remove(k);
		//		}

		@Override
		public String toString() {
			return additionOrder.toString();
		}

	}

	public static class NamedNode<T> extends TreeNode<T>{
		public String name;
		public NamedNode(String name, T data) {
			super(data);
			this.name = name;
		}

		public String toString() {
			return name+":"+data;
		}
	}

	public static class TreeNode<T> {

		public T data;
		public TreeNode<T> parent;
		public List<TreeNode<T>> children;

		public TreeNode() {
			super();
			children = new ArrayList<TreeNode<T>>();
		}

		public TreeNode(T data) {
			this();
			setData(data);
		}

		public List<TreeNode<T>> getChildren() {
			return this.children;
		}

		public int getNumberOfChildren() {
			return getChildren().size();
		}

		public boolean hasParent() {
			return parent != null;
		}

		public boolean hasChildren() {
			return (getNumberOfChildren() > 0);
		}

		public void setChildren(List<TreeNode<T>> children) {
			for(TreeNode<T> child : children){
				child.setParent(this);
			}
			this.children = children;
		}

		public void addChild(TreeNode<T> child) {
			child.setParent(this);
			children.add(child);
		}

		public void addChildAt(int index, TreeNode<T> child) throws IndexOutOfBoundsException {
			child.setParent(this);
			children.add(index, child);
		}

		public void setParent(TreeNode<T> parent){
			this.parent = parent;
		}

		public void removeChildren() {
			for(TreeNode<T> child : children){
				child.parent = null;
			}
			this.children = new ArrayList<TreeNode<T>>();
		}

		public void removeChildAt(int index) throws IndexOutOfBoundsException {
			children.get(index).parent=null;
			children.remove(index);
		}

		public TreeNode<T> getChildAt(int index) throws IndexOutOfBoundsException {
			return children.get(index);
		}

		public T getData() {
			return this.data;
		}

		public void setData(T data) {
			this.data = data;
		}

		public String toString() {
			return getData().toString();
		}

		public boolean equals(TreeNode<T> node) {
			return node.getData().equals(getData());
		}

		public int hashCode() {
			return getData().hashCode();
		}

		public String toStringVerbose() {
			String stringRepresentation = getData().toString() + ":[";

			for (TreeNode<T> node : getChildren()) {
				stringRepresentation += node.getData().toString() + ", ";
			}

			//Pattern.DOTALL causes ^ and $ to match. Otherwise it won't. It's retarded.
			Pattern pattern = Pattern.compile(", $", Pattern.DOTALL);
			Matcher matcher = pattern.matcher(stringRepresentation);

			stringRepresentation = matcher.replaceFirst("");
			stringRepresentation += "]";

			return stringRepresentation;
		}

		public String toNewickString() {
			String newickString = "";
			if(this.hasChildren()){
				newickString += "(";
				List<String> childNewickStrings = new LinkedList<String>();
				for(TreeNode<T> child : children){
					childNewickStrings.add(child.toNewickString());
				}
				newickString += StringUtils.join(childNewickStrings,",") + ")";
			}
			newickString += toString();
			return newickString;
		}
	}

	public static enum TreeTraversalOrderEnum { PRE_ORDER, POST_ORDER }

	public static class Tree<T> {

		private TreeNode<T> root;

		public Tree() {
		}

		public TreeNode<T> getRoot() {
			return this.root;
		}

		public void setRoot(TreeNode<T> root) {
			this.root = root;
		}

		public int getNumberOfNodes() {
			int numberOfNodes = 0;

			if(root != null) {
				numberOfNodes = auxiliaryGetNumberOfNodes(root) + 1; //1 for the root!
			}

			return numberOfNodes;
		}

		private int auxiliaryGetNumberOfNodes(TreeNode<T> node) {
			int numberOfNodes = node.getNumberOfChildren();

			for(TreeNode<T> child : node.getChildren()) {
				numberOfNodes += auxiliaryGetNumberOfNodes(child);
			}

			return numberOfNodes;
		}

		public boolean exists(TreeNode<T> nodeToFind) {
			return (find(nodeToFind) != null);
		}

		public TreeNode<T> find(TreeNode<T> nodeToFind) {
			TreeNode<T> returnNode = null;

			if(root != null) {
				returnNode = auxiliaryFind(root, nodeToFind);
			}

			return returnNode;
		}

		private TreeNode<T> auxiliaryFind(TreeNode<T> currentNode, TreeNode<T> nodeToFind) {
			TreeNode<T> returnNode = null;
			int i = 0;

			if (currentNode.equals(nodeToFind)) {
				returnNode = currentNode;
			}

			else if(currentNode.hasChildren()) {
				i = 0;
				while(returnNode == null && i < currentNode.getNumberOfChildren()) {
					returnNode = auxiliaryFind(currentNode.getChildAt(i), nodeToFind);
					i++;
				}
			}

			return returnNode;
		}

		public boolean isEmpty() {
			return (root == null);
		}

		public List<TreeNode<T>> build(TreeTraversalOrderEnum traversalOrder) {
			List<TreeNode<T>> returnList = null;

			if(root != null) {
				returnList = build(root, traversalOrder);
			}

			return returnList;
		}

		public List<TreeNode<T>> build(TreeNode<T> node, TreeTraversalOrderEnum traversalOrder) {
			List<TreeNode<T>> traversalResult = new ArrayList<TreeNode<T>>();

			if(traversalOrder == TreeTraversalOrderEnum.PRE_ORDER) {
				buildPreOrder(node, traversalResult);
			}

			else if(traversalOrder == TreeTraversalOrderEnum.POST_ORDER) {
				buildPostOrder(node, traversalResult);
			}

			return traversalResult;
		}

		private void buildPreOrder(TreeNode<T> node, List<TreeNode<T>> traversalResult) {
			traversalResult.add(node);

			for(TreeNode<T> child : node.getChildren()) {
				buildPreOrder(child, traversalResult);
			}
		}

		private void buildPostOrder(TreeNode<T> node, List<TreeNode<T>> traversalResult) {
			for(TreeNode<T> child : node.getChildren()) {
				buildPostOrder(child, traversalResult);
			}

			traversalResult.add(node);
		}

		public Map<TreeNode<T>, Integer> buildWithDepth(TreeTraversalOrderEnum traversalOrder) {
			Map<TreeNode<T>, Integer> returnMap = null;

			if(root != null) {
				returnMap = buildWithDepth(root, traversalOrder);
			}

			return returnMap;
		}

		public Map<TreeNode<T>, Integer> buildWithDepth(TreeNode<T> node, TreeTraversalOrderEnum traversalOrder) {
			Map<TreeNode<T>, Integer> traversalResult = new LinkedHashMap<TreeNode<T>, Integer>();

			if(traversalOrder == TreeTraversalOrderEnum.PRE_ORDER) {
				buildPreOrderWithDepth(node, traversalResult, 0);
			}

			else if(traversalOrder == TreeTraversalOrderEnum.POST_ORDER) {
				buildPostOrderWithDepth(node, traversalResult, 0);
			}

			return traversalResult;
		}

		private void buildPreOrderWithDepth(TreeNode<T> node, Map<TreeNode<T>, Integer> traversalResult, int depth) {
			traversalResult.put(node, depth);

			for(TreeNode<T> child : node.getChildren()) {
				buildPreOrderWithDepth(child, traversalResult, depth + 1);
			}
		}

		private void buildPostOrderWithDepth(TreeNode<T> node, Map<TreeNode<T>, Integer> traversalResult, int depth) {
			for(TreeNode<T> child : node.getChildren()) {
				buildPostOrderWithDepth(child, traversalResult, depth + 1);
			}

			traversalResult.put(node, depth);
		}

		public String toString() {
			/*
	        We're going to assume a pre-order traversal by default
			 */

			String stringRepresentation = "";

			if(root != null) {
				stringRepresentation = build(TreeTraversalOrderEnum.PRE_ORDER).toString();

			}

			return stringRepresentation;
		}

		public String toStringWithDepth() {
			/*
	        We're going to assume a pre-order traversal by default
			 */

			String stringRepresentation = "";

			if(root != null) {
				stringRepresentation = buildWithDepth(TreeTraversalOrderEnum.PRE_ORDER).toString();
			}

			return stringRepresentation;
		}

		public String toNewickString() {
			/*
	        We're going to assume a pre-order traversal by default
			 */

			String stringRepresentation = "(";

			if(root != null) {
				stringRepresentation += root.toNewickString();
			}

			stringRepresentation +=  ");";

			return stringRepresentation;
		}
	}

	public static <T> Set<T> difference(Set<T> set1, Set<T> set2){
		Set<T> union = Util.set(set1);
		union.addAll(set2);
		Set<T> intersection = intersection(set1,set2);
		union.removeAll(intersection);
		return union;
	}

	public static <T> Set<T> union(Set<T> set1, Set<T> set2){
		Set<T> union = Util.set(set1);
		union.addAll(set2);
		return union;
	}

	public static <T> Set<T> intersection(Set<T> set1, Set<T> set2){
		Set<T> intersection = Util.set(set1);
		intersection.retainAll(set2);
		return intersection;
	}

	public static <T> Set<T> subtract(Set<T> set1, Set<T> set2){
		Set<T> subtraction = Util.set(set1);
		subtraction.removeAll(set2);
		return subtraction;
	}

	public static <T extends Comparable> List<Integer> sortedOrder(List<T> c, boolean ascending){
		List<Integer> sorted = sort(seq(c.size()), c, true);
		List<Integer> order = new ArrayList<Integer>();
		for(int i=0; i<c.size(); i++){
			order.add(sorted.indexOf(i));
		}
		return order;
	}

	public static <T> List<T> enlist(T... objects){
		List<T> list = new ArrayList<T>();
		for(T t : objects){
			list.add(t);
		}
		return list;
	}

	public static <T> Map<T,Double> average( Map<T,Collection<? extends Number>> toAverage){
		Map<T,Double> average = new HashMap<T, Double>();
		for(Entry<T,Collection<? extends Number>> e : toAverage.entrySet()){
			average.put(e.getKey(), Util.mean(e.getValue()));
		}
		return average;
	}

	public static List<Double> rescale(List<? extends Number> numbers){
		Double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
		for(Number n : numbers){
			min = Math.min(n.doubleValue(), min);
			max = Math.max(n.doubleValue(), max);
		}
		double range = max-min;
		List<Double> rescaled = new ArrayList<Double>(numbers.size());
		for(Number n : numbers){
			rescaled.add((n.doubleValue()-min)/range);
		}
		return rescaled;
	}

	public static class countingIterator<E> implements Iterator<E>{
		Iterator<E> it;
		int i;
		public countingIterator(Iterator<E> it) {
			this.it = it;
			i=-1;
		}

		public int count(){
			return i;
		}

		public boolean hasNext() {
			return it.hasNext();
		}

		public E next() {
			i++;
			return it.next();
		}

		public void remove() {
			it.remove();
		}
	}

	public static void save(Object obj, String fileName) throws FileNotFoundException, IOException{
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(fileName));
		oos.writeObject(obj);
		oos.close();
	}
	
	public static void save(Object obj, OutputStream os) throws FileNotFoundException, IOException{
		ObjectOutputStream oos = new ObjectOutputStream(os);
		oos.writeObject(obj);
		oos.close();
	}

	public static Object load(String fileName) throws ClassNotFoundException, FileNotFoundException, IOException{
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fileName));
		return ois.readObject();
	}

	public static double logSum(List<Double> logNumbers){
		double max = Collections.max(logNumbers);

		double sum = 0;

		for(int i=0; i<logNumbers.size(); i++){
			double cur = logNumbers.get(i);
			if(Double.NEGATIVE_INFINITY == cur)
				continue;
			sum += Math.exp(cur - max);
		}

		return max + Math.log(sum);
	}
	
	/**
	 * 
	 * @param log_x
	 * @param log_y
	 * @return log(x-y)
	 */
	public static double logSubtract(double log_x, double log_y){
		double dif = Math.exp(log_x-log_y)-1;
		return Math.log(dif) + log_y;
	}

	public static double sum(Collection<? extends Number> numbers){
		double sum = 0;		

		for(Number n : numbers){
			sum += n.doubleValue();
		}

		return sum;
	}

	public static <T> void remove(Collection<T> collection, Collection<T> toRemove){
		for(T t : toRemove){
			collection.remove(t);
		}
	}

	public static <T> List<T> removeAll(Collection<T> collection, T toRemove){
		List<T> purged = Util.list(collection);
		while(purged.contains(toRemove))
			purged.remove(toRemove);
		return purged;
	}	

	public static double sse(Collection<? extends Number> numbers) {
		double sse = 0;
		double mean = mean(numbers);
		for(Number n : numbers){
			sse += Math.pow(mean-n.doubleValue(),2);
		}
		return sse;
	}
	
	public static double sum_squares(Collection<? extends Number> numbers) {
		double sum_squares = 0;
		for(Number n : numbers){
			double d = n.doubleValue();
			sum_squares += d*d;
		}
		return sum_squares;
	}

	public static double var(Collection<? extends Number> numbers) {
		double var = 0;
		double mean = mean(numbers);
		for(Number n : numbers){
			var += Math.pow(mean-n.doubleValue(),2);
		}
		return var/numbers.size();
	}

	public static double cov(List<? extends Number> numbers_x, List<? extends Number> numbers_y) {
		double cov = 0;
		double mean_x = mean(numbers_x);
		double mean_y = mean(numbers_y);
		for(int i=0; i<numbers_x.size(); i++){
			cov += (numbers_x.get(i).doubleValue()-mean_x)*(numbers_y.get(i).doubleValue()-mean_y);
		}
		cov /= numbers_x.size();
		return cov;
	}
	
	/**
	 * entropy base e 
	 * @param p_i
	 * @return
	 */
	public static double entropy(List<Double> p_i){
		double entropy = 0;
		
		for(Double d : p_i){
			if(d > 0)
				entropy += d*Math.log(d);
		}
		
		return -entropy;
	}
	
	/**
	 * For entropy base e the value of JSD will be between 0 and ln(2)
	 * @param p_i, a set of discrete distributions with a common domain
	 * @return the Jensen Shannon divergence between the distributions
	 */
	public static double JSD(List<Double>... p_i){
		int l = p_i[0].size();
		int n = p_i.length;
		double pi = 1./n;
		List<Double> p_mu = new ArrayList<Double>(l);
		
		// JSD is equal to the difference between the entropy
		// of the average distribution, and the average entropy
		// of the individual distributions
		for(int i=0; i<l; i++){
			double p_mu_i = 0;
			for(int j=0; j<n; j++){
				p_mu_i += pi*p_i[j].get(i);
			}
			p_mu.add(p_mu_i);
		}
		
		double JSD = entropy(p_mu);
		for(int j=0; j<n; j++){
			JSD -= pi*entropy(p_i[j]);
		}
		
		return JSD;
	}

	public static <S extends Number, T extends Number> double pearson_correlation(List<S> numbers_x, List<T> numbers_y) {
		double cor = cov(numbers_x,numbers_y);
		double stdev_x = stdev(numbers_x);
		double stdev_y = stdev(numbers_y);
		cor /= stdev_x*stdev_y;
		return cor;
	}
	
	public static <T extends Comparable<T>> List<Double> rank(List<T> to_rank){
		List<Integer> order = sort(Util.seq(to_rank.size()), to_rank, true);
		
		List<Double> ranks = new ArrayList<Double>(to_rank.size());

		int n=1;
		T last=to_rank.get((int) order.get(0));
		double sum_rank_y=0;

		for(int i=1; i<to_rank.size(); i++){
			T next = to_rank.get((int) order.get(i));
			
			if(next.compareTo(last)==0){
				sum_rank_y+=i;
				n++;
			}
			else{
				double avg_rank = sum_rank_y/n;
				for(int j=0; j<n; j++){
					ranks.add(avg_rank);
				}
				sum_rank_y=i;
				n=1;
				last=next;
			}
		}
		
		if(n!=0){
			double avg_rank = sum_rank_y/n;
			for(int j=0; j<n; j++){
				ranks.add(avg_rank);
			}
		}
		
	List<Double> y_rank = repeat(null, to_rank.size());

//	System.out.println(ranks);
		for(int i=0; i<order.size(); i++){
			y_rank.set(order.get(i), ranks.get(i));
		}
//		System.out.println(y_rank);
		return y_rank;
	}

	public static <S extends Comparable<S>, T extends Comparable<T>> double spearman_correlation(List<S> numbers_x, List<T> numbers_y) {
		List<Integer> x_order = sort(Util.seq(numbers_x.size()), numbers_x, true);
		List<Integer> y_order = sort(Util.seq(numbers_y.size()), numbers_y, true);

		List<S> x_sorted = repeat(null, numbers_x.size());
		List<T> y_sorted = repeat(null, numbers_y.size());
	
		for(int i=0; i<x_order.size(); i++){
			x_sorted.set(i, numbers_x.get((int) x_order.get(i)));
			y_sorted.set(i, numbers_y.get((int) y_order.get(i)));
		}

		List<Double> x_ranks = new ArrayList<Double>(numbers_x.size());
		List<Double> y_ranks = new ArrayList<Double>(numbers_x.size());

		int n_x=1;
		int n_y=1;
		S last_x=x_sorted.get(0);
		T last_y=y_sorted.get(0);
		double sum_rank_x=0;
		double sum_rank_y=0;

		for(int i=1; i<x_sorted.size(); i++){
			S next_x = x_sorted.get(i);
			T next_y = y_sorted.get(i);

//			System.out.printf("next %s\n", next_x);
//			System.out.println(next_x.compareTo(last_x)==0);

			if(next_x.compareTo(last_x)==0){
				sum_rank_x+=i;
				n_x++;
			}
			else{
				double avg_rank = sum_rank_x/n_x;
				for(int j=0; j<n_x; j++){
					x_ranks.add(avg_rank);
				}
				sum_rank_x=i;
				n_x=1;
				last_x=next_x;
			}
			if(next_y.compareTo(last_y)==0){
				sum_rank_y+=i;
				n_y++;
			}
			else{
				double avg_rank = sum_rank_y/n_y;
				for(int j=0; j<n_y; j++){
					y_ranks.add(avg_rank);
				}
				sum_rank_y=i;
				n_y=1;
				last_y=next_y;
			}
		}
		if(n_x!=0){
			double avg_rank = sum_rank_x/n_x;
			for(int j=0; j<n_x; j++){
				x_ranks.add(avg_rank);
			}
		}
		if(n_y!=0){
			double avg_rank = sum_rank_y/n_y;
			for(int j=0; j<n_y; j++){
				y_ranks.add(avg_rank);
			}
		}
//		System.out.println(x_sorted);
//		System.out.println(y_sorted);
//		
//		System.out.println(x_order);
//		System.out.println(y_order);
//
//		System.out.println(x_ranks);
//		System.out.println(y_ranks);
		
		List<Double> x_rank = repeat(null, numbers_x.size());
		List<Double> y_rank = repeat(null, numbers_y.size());

		for(int i=0; i<x_order.size(); i++){
			x_rank.set(x_order.get(i), x_ranks.get(i));
			y_rank.set(y_order.get(i), y_ranks.get(i));
		}

//		System.out.println(x_rank);
//		System.out.println(y_rank);
		
		return pearson_correlation(x_rank, y_rank);
	}

	public static double mutual_information(List<Double> prob_xy, List<Double> prob_x, List<Double> prob_y) {

		return 0;
	}

	public static double r_squared(List<? extends Number> numbers_x, List<? extends Number> numbers_y) {
		double cor = pearson_correlation(numbers_x,numbers_y);	
		return cor*cor;
	}

	public static double stdev(Collection<? extends Number> numbers) {
		return Math.sqrt(var(numbers));
	}

	public static double stderr(Collection<? extends Number> numbers) {
		return Math.sqrt(var(numbers)/numbers.size());
	}

	public static double mean(Collection<? extends Number> numbers){
		return sum(numbers)/numbers.size();
	}

	/**
	 * 
	 * @param numbers collection of numbers
	 * @param quantile
	 * @return
	 */
	public static <T extends Comparable<? super T>> T quantile(Collection<T> numbers, double quantile){
		List<T> sorted = new ArrayList<T>();
		for(T n : numbers){
			sorted.add(n);
		}
		Collections.sort(sorted);

		int index = (int) Math.ceil(quantile*(numbers.size()-1));
		return sorted.get(index);
	}

	public static double median(Collection<? extends Number> numbers){
		List<Double> sorted = new ArrayList<Double>();
		for(Number n : numbers){
			sorted.add(n.doubleValue());
		}
		Collections.sort(sorted);
		int size = sorted.size();
		if(size % 2 == 1)
			return sorted.get(size/2);
		else
			return (sorted.get((size-1)/2)+sorted.get(size/2))/2;
	}

	public static double product(Collection<? extends Number> numbers){
		double product = 1;
		Number[] nums = numbers.toArray(new Number[0]);

		for(Number n : nums){
			product *= n.doubleValue();
		}
		return product;
	}

	public static List<Byte> list(byte[] bytes){
		return list(ArrayUtils.toObject(bytes));
	}

	public static List<Short> list(short[] shorts){
		return list(ArrayUtils.toObject(shorts));
	}

	public static List<Integer> list(int[] ints){
		return list(ArrayUtils.toObject(ints));
	}

	public static List<Long> list(long[] longs){
		return list(ArrayUtils.toObject(longs));
	}

	public static List<Float> list(float[] floats){
		return list(ArrayUtils.toObject(floats));
	}

	public static List<Double> list(double[] doubles){
		return list(ArrayUtils.toObject(doubles));
	}

	public static List<Character> list(char[] chars){
		return list(ArrayUtils.toObject(chars));
	}

	public static List<Boolean> list(boolean[] booleans){
		return list(ArrayUtils.toObject(booleans));
	}
	
	public static <T> Queue<T> queue(T... objects){
		return new LinkedList<T>(Util.list(objects));
	}

	public static <T> List<T> list(T... objects){
		List list = new ArrayList(objects.length);
		for(T t : objects){
			if(t != null && t.getClass().isArray()){
				if(!t.getClass().getComponentType().isArray()){
					List primitiveList = null;
					if(objects.length == 1)
						primitiveList = list;
					else
						primitiveList = new ArrayList(Array.getLength(t));
					for(int i=0;i<Array.getLength(t);i++)
						primitiveList.add(Array.get(t, i));
					if(objects.length != 1)
						list.add(primitiveList);					
				}
				else
					list.add(list((T[]) t));
			}
			else
				list.add(t);
		}
		return list;
	}

	public static <T> List<T> list(Collection<T> elements){
		ArrayList<T> list = new ArrayList<T>(elements);
		return list;
	}

	public static <T> Set<T> set(T ...elements){
		Set<T> set = new HashSet<T>(Arrays.asList(elements));
		return set;
	}

	public static <T> Set<T> set(Collection<T> elements){
		Set<T> set = new HashSet<T>(elements);
		return set;
	}

	public static List list(){
		return new ArrayList();
	}

	public static <T> int baseValue(List<T> values, List<T> bitString){
		int value = 0;
		int base = values.size();
		for(T bit : bitString){
			int i = values.indexOf(bit);
			value *= base;
			value += i;
		}
		return value;
	}

	/**
	 * Makes no assumptions about the map
	 * @param map
	 * @return
	 */
	public static <S,T> Map<S,Collection<T>> mapInverse(Map<T,S> map){
		Map<S,Collection<T>> inverse = new HashMap<S, Collection<T>>();
		for(Entry<T, S> e : map.entrySet()){
			if(!inverse.containsKey(e.getValue()))
				inverse.put(e.getValue(), new LinkedList<T>());
			inverse.get(e.getValue()).add(e.getKey());
		}
		return inverse;
	}

	/**
	 * Assumes that the mappings are 1:1, any n:1 mappings will be lost
	 * @param map
	 * @return
	 */
	public static <S,T> Map<S,T> mapInvert(Map<T,S> map){
		Map<S,T> inverse = new HashMap<S,T>();
		for(Entry<T, S> e : map.entrySet()){

			inverse.put(e.getValue(),e.getKey());
		}
		return inverse;
	}

	public static <T> List<T> select(T[] array,	List<Integer> indices) {
		ArrayList<T> selected = new ArrayList<T>(indices.size());
		for(Integer i : indices){
			selected.add(array[i]);
		}
		return selected;
	}

	public static <T> int occurencesOf(Collection<T> c, T obj) {
		int occurences = 0;
		for(T t : c){
			if(t.equals(obj))
				occurences++;
		}
		return occurences;
	}

	public static <T> List<Integer> indicesOf(List<T> list, T obj) {
		List<Integer> indices = new LinkedList<Integer>();
		for(Integer i = 0 ; i < list.size() ; i++){
			if(list.get(i).equals(obj))
				indices.add(i);
		}
		return indices;
	}

	public static <T> List<T> select(List<T> list,	List<Integer> indices) {
		LinkedList<T> selected = new LinkedList<T>();
		for(Integer i : indices){
			selected.add(list.get(i));
		}
		return selected;
	}

	public static <T> Collection<T> combine(Collection<T> ...collections){
		Collection<T> combined = new ArrayList<T>();
		for(Collection<T> collection : collections){
			combined.addAll(collection);
		}
		return combined;
	}

	public static <T> Collection<T> combine(Collection<? extends Collection<T>> collections){
		Collection<T> combined = new ArrayList<T>();
		for(Collection<T> collection : collections){
			combined.addAll(collection);
		}
		return combined;
	}

	public static <T> List<T> reverse(List<T> toReverse){
		List<T> reversed = new ArrayList(toReverse);
		Collections.reverse(reversed);
		return reversed;
	}

	public static <T extends Collection<S>, S> T copy(T t){
		try {
			T copy = (T) t.getClass().newInstance();
			for(S s : t){
				copy.add(s);
			}
			return copy;
		} catch (Exception e) {
			e.printStackTrace();
		} 
		throw new RuntimeException("Error copying collection");
	}

	public static <T> List<T> cast(Class<? extends T> clazz, Collection<?> c) {
		List<T> r = new ArrayList<T>(c.size());
		for(Object o: c){
			r.add(clazz.cast(o));
		}
		return r;
	}

	public static <S,T> Map<S, T> map(List<S> keys, List<T> values) {
		Map<S,T> m = new HashMap<S, T>();
		for(int i=0; i<keys.size(); i++)
			m.put(keys.get(i), values.get(i));

		return m;
	}

	public static <T extends Object> List<List<T>> enumerateObjects(List<T> possibleValues, int listSize){

		List<List<T>> types = new LinkedList<List<T>>();

		for(T o : possibleValues){
			if(listSize > 1)
				for(List<T> remaining : enumerateObjects(possibleValues, listSize-1)){
					LinkedList<T> type = new LinkedList<T>();
					type.add(o);
					type.addAll(remaining);
					types.add(type);
				}
			else{
				LinkedList<T> type = new LinkedList<T>();
				type.add(o);
				types.add(type);
			}
		}

		return types;
	}


	public static class PartialIterator<T> implements Iterable<T>{

		int n;
		Iterable<T> iterable;

		public PartialIterator(Iterable<T> toIterate) {
			n=-1;
			iterable = toIterate;
		}

		public PartialIterator(Iterable<T> toIterate, int n) {
			this.n = n;
			iterable = toIterate;
		}

		class internalIterator implements Iterator<T>{
			int i,n;
			Iterator<T> iterator;

			public internalIterator(Iterable<T> toIterate, int n) {
				i=0;
				this.n = n;
				iterator = toIterate.iterator();
			}

			public boolean hasNext() {
				return iterator.hasNext() && (n == -1 || i < n);
			}

			public T next() {
				if(n != -1 && i>n-1)
					throw new NoSuchElementException();
				i++;
				return iterator.next();
			}

			public void remove() {
				iterator.remove();
			}
		}

		public Iterator<T> iterator() {
			return new internalIterator(iterable, n);
		}

	}

	public static class Count implements Iterable<Integer>, Iterator<Integer>{

		int size;
		int i;

		public Count(int size){
			this.size = size;	
			i=0;
		}


		public Iterator<Integer> iterator() {			
			return this;
		}


		public boolean hasNext() {
			return i<size;
		}


		public Integer next(){
			return i++;
		}

		public void remove() {	
		}

	}

	public static class Timer{
		long startTime;
		public Timer(){
			restart();
		}

		public void restart() {
			startTime = System.nanoTime();
		}

		public long elapsedNanos(){
			return System.nanoTime()-startTime;
		}

		public double elapsedMillis(){
			return elapsedNanos()/1000000.0;
		}

		public double elapsedTime(){
			return elapsedNanos()/1000000000.0;
		}

	}

	public static class multiIterator<E> implements Iterator<E[]>, Iterable<E[]>{
		Iterator<E>[] iterators;
		boolean nullOnEmpty;
		Class<E> type;

		@SuppressWarnings("unchecked")
		public multiIterator(boolean nullOnEmpty, Class<E> type, Iterable<E> ...iterables){
			this.type = type;
			this.nullOnEmpty=nullOnEmpty;
			iterators = new Iterator[iterables.length];
			for (int i = 0; i < iterables.length; i++) {
				iterators[i] = iterables[i].iterator();
			}		
		}


		public boolean hasNext() {
			boolean hasNext = false;
			for(Iterator<E> it : iterators){
				hasNext |= it.hasNext();
			}
			return hasNext;
		}


		@SuppressWarnings("unchecked")
		public E[] next() {
			E[] next = (E[]) Array.newInstance(type, iterators.length);
			for (int i = 0; i < iterators.length; i++) {
				if(nullOnEmpty && !iterators[i].hasNext())
					continue;
				else
					next[i] = iterators[i].next();
			}

			return next;
		}


		public void remove() {

		}


		public Iterator<E[]> iterator() {
			return this;
		}
	}

	public static class MultiIterable implements Iterable<Object[]>{
		Iterable[] iterables;
		boolean nullOnEmpty;
		int l;

		@SuppressWarnings("unchecked")
		public MultiIterable(boolean nullOnEmpty, Iterable ...iterables){
			this.nullOnEmpty=nullOnEmpty;
			l = iterables.length;
			this.iterables = iterables;
		}

		public Iterator<Object[]> iterator() {
			return new MultiIterator();
		}

		class MultiIterator implements Iterator<Object[]>{

			Iterator[] iterators;

			public MultiIterator(){

				iterators = new Iterator[l];
				for (int i = 0; i < l; i++) {
					iterators[i] = iterables[i].iterator();
				}		
			}


			public boolean hasNext() {
				boolean hasNext = false;
				for(Iterator it : iterators){
					hasNext |= it.hasNext();
				}
				return hasNext;
			}


			@SuppressWarnings("unchecked")
			public Object[] next() {
				Object[] next = new Object[l];
				for (int i = 0; i < l; i++) {
					if(nullOnEmpty && !iterators[i].hasNext())
						continue;
					else
						next[i] = iterators[i].next();
				}

				return next;
			}


			public void remove() {

			}
		}
	}

	public static class SequentialIterator<E> implements Iterable<E>{
		Iterable<E>[] iterables;

		public SequentialIterator(Iterable<E> ...iterables){
			this.iterables = new Iterable[iterables.length];
			for (int j = 0; j < iterables.length; j++) {
				this.iterables[j] = iterables[j];
			}		
		}

		public <T extends Iterable<E>> SequentialIterator(Collection<T> iterables){
			this.iterables = new Iterable[iterables.size()];
			Iterator<T> iterablesIterator = iterables.iterator();
			for (int j = 0; j < iterables.size(); j++) {
				this.iterables[j] = iterablesIterator.next();
			}		
		}

		public Iterator<E> iterator() {
			return new IterableIterator(iterables);
		}

		private class IterableIterator implements Iterator<E> {
			Iterator<E>[] iterators;
			int i;

			public IterableIterator(Iterable<E> ...iterables){
				iterators = new Iterator[iterables.length];
				i=0;
				for (int j = 0; j < iterables.length; j++) {
					iterators[j] = iterables[j].iterator();
				}		
				increment();
			}

			private void increment(){
				while(i < iterators.length && !iterators[i].hasNext()){
					i++;
				}
			}


			public boolean hasNext() {
				return i < iterators.length && iterators[i].hasNext();
			}


			public E next() {
				E next = iterators[i].next();
				increment();			
				return next;
			}

			public void remove() {

			}
		}
	}

	public static class  comboIterator<T> implements Iterator<List<T>>, Iterable<List<T>> {
		List<T> values;
		int numValues;
		int[] iterators;
		int numElements;
		List<T> last;
		int numCombos;
		int comboCount;
		boolean hasNext;

		/**
		 * 
		 * @param possibleValues
		 * @param listSize - must be greater than 0
		 */
		public comboIterator(List<T> possibleValues, int listSize) {
			values = possibleValues;
			numValues = possibleValues.size();

			iterators = new int[listSize];
			Arrays.fill(iterators, 1);
			iterators[listSize-1] = 0;

			numElements = listSize;

			last = new ArrayList<T>(listSize);
			for(int i=0; i<numElements; i++){
				last.add(possibleValues.get(0));
			}
			numCombos = (int) Math.pow(numValues, numElements);
			hasNext = true;
		}


		public boolean hasNext() {
			return hasNext;
		}


		public List<T> next() {
			if(!hasNext())
				throw new NoSuchElementException();
			for(int currentIndex=numElements-1; currentIndex>-1; currentIndex--){
				if(iterators[currentIndex] < numValues){
					last.set(currentIndex,values.get(iterators[currentIndex]++));
					break;
				}
				else{
					iterators[currentIndex] = 0;
					last.set(currentIndex, values.get(iterators[currentIndex]++));
				}
			}
			comboCount++;
			hasNext = comboCount < numCombos;
			return last;
		}


		public void remove() {}


		public Iterator<List<T>> iterator() {
			return this;
		}
	};

	public static class ComboIndexer<T>{
		List<T> values;
		int radix;

		public ComboIndexer(List<T> values) {
			if(values.size() != Util.set(values).size())
				throw new IllegalArgumentException("All values must be distinct");
			this.values = values;
			radix = values.size();
		}

		public long indexOf(List<T> toIndex){
			long index = 0;

			for(T t : toIndex){
				index *= radix;
				index += values.indexOf(t);
			}

			return index;
		}
	}

	/**
	 * Iterates over all possible subsets of a set of values. There are n choose x different subsets.
	 *
	 * @param <T> - The type of the value that the collection is made of 
	 */
	public static class  subsetIterator<T> implements Iterator<Collection<T>>, Iterable<Collection<T>> {
		List<T> values;
		int numValues;
		int[] iterators;
		int setSize;
		List<T> last;
		long numCombos;
		int comboCount;
		boolean hasNext;

		/**
		 * 
		 * @param possibleValues
		 * @param subsetSize - must be greater than 0
		 */
		public subsetIterator(Collection<T> possibleValues, int subsetSize) {
			values = Util.list(possibleValues);
			numValues = possibleValues.size();

			iterators = new int[subsetSize];

			setSize = subsetSize;

			last = new ArrayList<T>(subsetSize);
			for(int i=0; i<setSize; i++){
				iterators[i] = i-1;
				last.add(values.get(i));
			}

			numCombos = ArithmeticUtils.binomialCoefficient(possibleValues.size(), subsetSize);
			hasNext = true;
		}


		public boolean hasNext() {
			return hasNext;
		}


		public List<T> next() {
			if(!hasNext())
				throw new NoSuchElementException();
			for(int currentIndex=0; currentIndex < setSize; currentIndex++){
				iterators[currentIndex]++;
				if(currentIndex+1 < setSize && iterators[currentIndex] == iterators[currentIndex+1]){
					iterators[currentIndex] = currentIndex;
					last.set(currentIndex, values.get(currentIndex));
				}
				else{
					last.set(currentIndex, values.get(iterators[currentIndex]));
					break;
				}			
			}
			comboCount++;
			hasNext = comboCount < numCombos;
			return last;
		}


		public void remove() {}


		public Iterator<Collection<T>> iterator() {
			return this;
		}
	};

	public static class  PermutationIterator<T> implements Iterator<List<T>>, Iterable<List<T>> {
		List<T> toPermute;
		private int[] currentPermutation;         //current perm
		private int permutationLength;              //length of perm string ("N")
		private long totalPermutations;
		private long counter = 0;

		public PermutationIterator(List<T> toPermute) {
			this.toPermute = toPermute;
			permutationLength = toPermute.size();

			totalPermutations = ArithmeticUtils.factorial(permutationLength);
			currentPermutation = new int[permutationLength];
			for (int i=0; i < permutationLength; i++) {
				currentPermutation[i] = i;
			}
		}

		public long totalPermutations() {
			return totalPermutations;
		}

		public boolean hasNext() {
			return counter < totalPermutations;
		}    

		/*
		 * Algorithm by E. W. Dijkstra, "A Discipline of Programming",
		 *       Prentice-Hall, 1976, p.71
		 * Another example (Public Domain) by Tim Tyler from 2000 is at:
		 *       http://mandala.co.uk/permutations/
		 *       http://mandala.co.uk/permutations/plugin/Permutations.java
		 */
		public List<T> next() {

			counter++;
			if (counter == 1) {
				return Util.select(toPermute, Util.list(currentPermutation));
			} else if (counter > totalPermutations) {
				throw new NoSuchElementException();
			} 

			int i = permutationLength - 1;

			while(currentPermutation[i - 1] >= currentPermutation[i]) {
				i--;
			}

			int j = permutationLength;

			while(currentPermutation[j - 1] <= currentPermutation[i - 1]) {
				j--;
			}

			swap(i - 1, j - 1);

			i++;
			j = permutationLength;

			while(i < j) {
				swap(i - 1, j - 1);
				i++;
				j--;
			}

			return Util.select(toPermute, Util.list(currentPermutation));
		}

		//abstract remove() specified by Iterator is not implemented
		public void remove() {}

		private void swap(int index1, int index2) {
			int tmp;
			tmp = currentPermutation[index1];
			currentPermutation[index1] = currentPermutation[index2];
			currentPermutation[index2] = tmp;
		}

		public Iterator<List<T>> iterator() {
			return new PermutationIterator<T>(toPermute);
		}
	}

	public static class BitString implements Iterable<Integer>{
		int size, arraySize, overflow;
		int[] bits;
		final int SIZE; 

		public BitString(int size) {
			SIZE = Integer.SIZE;
			this.size = size;
			overflow = size%SIZE;
			arraySize = size/SIZE + (overflow >0 ? 1 : 0);
			bits = new int[arraySize];
		}

		public void set(int bit){
			bits[bit/SIZE] |= 1 << (bit % SIZE); 
		}

		public void unset(int bit){
			bits[bit/SIZE] &= ~(1 << (bit % SIZE)); 
		}

		public int get(int bit){
			return 1 & (bits[bit/SIZE] >> bit); 
		}

		public BitString and(BitString other){
			BitString maxBits = arraySize > other.arraySize ? this : other;
			BitString minBits = maxBits == this ? other : this;
			int minSize = minBits.arraySize;
			BitString and = new BitString(maxBits.size);
			for(int i=0; i<minSize; i++){
				and.bits[i] = maxBits.bits[i] & minBits.bits[i];
			}
			return and;		
		}

		public BitString or(BitString other){
			BitString maxBits = arraySize > other.arraySize ? this : other;
			BitString minBits = maxBits == this ? other : this;
			int maxSize = maxBits.arraySize;
			int minSize = minBits.arraySize;
			BitString or = new BitString(maxBits.size);
			for(int i=0; i<minSize; i++){
				or.bits[i] = maxBits.bits[i] | minBits.bits[i];
			}
			for(int i=minSize; i<maxSize; i++){
				or.bits[i] = maxBits.bits[i];
			}
			return or;		
		}
		
		public void OR(BitString other){
			for(int i=0; i<arraySize; i++){
				bits[i] = bits[i] | other.bits[i];
			}		
		}

		public BitString not(){
			BitString not = new BitString(size);
			int mask = ~(~0 << overflow);
			for(int i=0; i<arraySize; i++){
				not.bits[i] = ~bits[i];
			}
			not.bits[arraySize-1] &= mask;
			return not;		
		}

		public BitString xor(BitString other){
			BitString maxBits = arraySize > other.arraySize ? this : other;
			BitString minBits = maxBits == this ? other : this;
			int maxSize = maxBits.arraySize;
			int minSize = minBits.arraySize;
			BitString xor = new BitString(maxBits.size);

			for(int i=0; i<minSize; i++){
				xor.bits[i] = maxBits.bits[i] ^ minBits.bits[i];
			}
			for(int i=minSize; i<maxSize; i++){
				xor.bits[i] = maxBits.bits[i] ^ 0;
			}

			return xor;	
		}

		public String toString(){
			StringBuffer sb = new StringBuffer(size);
			for(int i=size-1; i>=0; i--){
				sb.append(get(i));
			}
			return sb.toString();
		}
		
		public int count(){
			int count = 0;
			for(int i=0; i<arraySize; i++){
				count += Integer.bitCount(bits[i]);
			}
			return count;
		}

		public Iterator<Integer> iterator() {
			final int s = size;
			final int[] b = bits; 
			return new Iterator<Integer>() {
				int i=-1, arrayI=-1;
				int cur;

				public boolean hasNext() {				
					return i < s-1;
				}

				public Integer next() {
					i++;
					if(i % SIZE == 0){
						arrayI++;
						cur = b[arrayI];
					}
					int next = cur & 1;
					cur >>= 1;
					return next;
				}

				public void remove() {
					//do nothing
				}
			};
		}
	}

	/**
	 * Samples the states using a uniform distribution
	 * @param states - The states to be sampled
	 * @return - The sampled state
	 */
	public static <T> T sample(Collection<T> states) {
		List<T> samples = new ArrayList<T>(states);
		return samples.get((int) (Math.random()*samples.size()));
	}

	public static <T> List<T> sample(List<Double> probabilities, List<T> states, int number){
		List<Double> cumulative = new ArrayList<Double>(probabilities.size());
		List<T> samples = new ArrayList<T>(number);
		double sum = 0;
		for(Double p : probabilities){
			sum += p;
			cumulative.add(sum);
		}
		for(int i=0; i<number; i++){
			int index = Collections.binarySearch(cumulative, sum*Math.random());
			if(index < 0)
				index = -index -1;

			samples.add(states.get(index));
		}

		return samples;		
	}

	public static <T> T sample(List<Double> probabilities, List<T> states){
		List<Double> cumulative = new ArrayList<Double>(probabilities.size());
		T sample = null;
		double sum = 0;
		List<T> toSample = new ArrayList<T>();
		for(int i=0; i< probabilities.size(); i++){
			Double p = probabilities.get(i);
			T t = states.get(i);
			sum += p;
			if(p>0){
				toSample.add(t);
				cumulative.add(sum);
			}
		}
		int index = Collections.binarySearch(cumulative, sum*Math.random());
		if(index < 0)
			index = -index -1;

		sample = (toSample.get(index));

		return sample;		
	}

	public static List<Double> normalize(List<? extends Number> probabilities){
		List<Double> normalized = new ArrayList<Double>();
		double sum = 0;
		for(Number n : probabilities){
			sum += n.doubleValue();
		}
		for(Number n : probabilities){
			normalized.add(n.doubleValue()/sum);
		}
		return normalized;		
	}
	
	public static List<Double> standardize(List<? extends Number> numbers){
		List<Double> standardized = new ArrayList<Double>(numbers.size());
		double mean = Util.mean(numbers);
		double stdev = Util.stdev(numbers);
		for(Number n : numbers){
			standardized.add((n.doubleValue()-mean)/stdev);
		}
		return standardized;		
	}

	public static List<Double> scale(List<? extends Number> numbers, double min, double max){
		List<Double> scaled = new ArrayList<Double>(numbers.size());
		double MIN = numbers.get(0).doubleValue();
		double MAX = MIN;

		for(Number n : numbers){
			double d = n.doubleValue();
			if(d > MAX)
				MAX = d;
			if(d < MIN)
				MIN = d;
		}
		if(MAX==MIN){
			throw new RuntimeException("max cannot equal min");
		}
		for(Number n : numbers){
			double d = n.doubleValue();
			scaled.add(((d-MIN)/(MAX-MIN)*(max-min)) + min);
		}
		return scaled;		
	}

	public static List<Integer> seq(int j) {
		ArrayList<Integer> seq = new ArrayList<Integer>(j);
		for(Integer i=0; i<j; i++){
			seq.add(i);
		}
		return seq;
	}
	
	public static List<Integer> seq(int start, int end) {
		ArrayList<Integer> seq = new ArrayList<Integer>(end-start+1);
		for(Integer i=start; i<end; i++){
			seq.add(i);
		}
		return seq;
	}

	public static Double floor(Double[] items, double key) {
		int location = Arrays.binarySearch(items, key);

		if (location >= 0) return items[location];

		location = -(location + 1);

		if (location == 0) {
			return items[0];
		} else if (location == items.length) {
			return items[items.length - 1];
		}

		return items[location - 1];
	}

	public static Double ceiling(Double[] items, double key) {
		int location = Arrays.binarySearch(items, key);

		if (location >= 0) return items[location];

		location = -(location + 1);

		if (location == 0) {
			return items[0];
		} else if (location == items.length) {
			return items[items.length - 1];
		}

		return items[location];
	}

	public static void printf(String string, Object... objects) {
		String[] labels = string.split(",");
		for(int i=0; i<labels.length; i++){
			System.out.printf("%s : %s\n", labels[i].trim(), objects[i]);
		}
	}

	public static <T> T getObj(Collection<T> collection, T toGet){
		for(T t : collection){
			if(t.equals(toGet))
				return t;
		}
		throw new NoSuchElementException();
	}
	public static <S,T>  List<T> get(Map<S, T> map, Collection<S> list) {
		List<T> gotten = new ArrayList<T>(list.size());

		for(S s : list){
			gotten.add(map.get(s));
		}

		return gotten;
	}

	public static double logit(double x, double center, double scale) {
		return 1/(1 + Math.exp(scale*(x - center)));
	}

	public static String sprintf(String format, Object... toPrint){
		return String.format(format, toPrint);
	}

	public static <T> T get(Collection<T> collection, T toGet){
		for(T t : collection){
			if(t.equals(toGet))
				return t;
		}
		throw new RuntimeException("Collection does not contain any element eqal to " + toGet);
	}

	public static <S extends Comparable<S>> List<Integer> rank(List<S> values, boolean ascending) {
		return sort(seq(values.size()),values,ascending);
	}

	public static <S extends Comparable<S>,T> List<T> sort(List<T> listToSort, List<S> values, boolean ascending) {
		List<Tuple<S,T>> s = sortTuple(listToSort, values, ascending);
		List<T> sorted = new ArrayList<T>(s.size());
		for(int i=0; i<s.size(); i++){
			sorted.add(s.get(i).t);
		}
//		List<T> sorted = new ArrayList<T>();
//		List<S> order = Util.list(values);
//		List<S> list = new LinkedList<S>(values);
//
//		Collections.sort(order);
//
//		Iterator<S> it = order.iterator();
//	
//		while(it.hasNext()){
//			S s = it.next();
//
//			int i = list.indexOf(s);
//
//			if(ascending)
//				sorted.add(listToSort.get(i));
//			else
//				sorted.add(0, listToSort.get(i));
//
//			list.set(i, null);
//		}

		return sorted;
	}
	
	public static class Tuple<S,T> {
		public S s;
		public T t;
		
		public Tuple(S s, T t) {
		this.s = s;
		this.t = t;
		}
	}
	
	public static <S extends Comparable<S>,T> List<Tuple<S,T>> sortTuple(List<T> listToSort, List<S> values, boolean ascending) {
		
		Comparator<Tuple<S,T>> c = ascending ? new Comparator<Tuple<S,T>>() {
			public int compare(Tuple<S,T> o1, Tuple<S,T> o2) {
				return o1.s.compareTo(o2.s);
			}
		} : new Comparator<Tuple<S,T>>() {
			public int compare(Tuple<S,T> o1, Tuple<S,T> o2) {
				return o2.s.compareTo(o1.s);
			}
		};
		
		int n = listToSort.size();
		
		List<Tuple<S,T>> toSort = new ArrayList<Tuple<S,T>>(n);
		for(int i=0; i<n; i++){
			toSort.add(new Tuple<S,T>(values.get(i), listToSort.get(i)));
		}
		
		Collections.sort(toSort, c);

		return toSort;
	}

	public static Color hsba(Number h, Number s, Number b, Number a) {
		Color c = new Color(Color.HSBtoRGB(h.floatValue(),s.floatValue(),b.floatValue()));
		float[] components = c.getColorComponents(null);
		return new Color(components[0],components[1],components[2],a.floatValue());
	}

	public static Color hsb(Number h, Number s, Number b) {
		return new Color(Color.HSBtoRGB(h.floatValue(),s.floatValue(),b.floatValue()));
	}

	public static boolean notNull(Object... objects) {
		boolean notNull = true;
		for(Object obj : objects)
			notNull &= obj!=null;
		return notNull;
	}

	public static <T> List<T> repeat(T t, int size) {
		List<T> repeated = new ArrayList<T>(size);
		for(int i=0; i<size; i++){
			repeated.add(t);
		}
		return repeated;
	}

	public static abstract class Factory<T> {
		public abstract T create();
	}

	public static class FactoryMap<S,T>{
		Factory<T> factory;
		Map<S,T> map;

		public FactoryMap(Factory<T> factory) {
			this.factory = factory;
			map = new HashMap<S,T>();
		}

		public T get(S key){
			if(!map.containsKey(key)){
				map.put(key, factory.create());
			}
			return map.get(key);
		}

		public Collection<S> keySet(){
			return map.keySet();
		}
	}

	public static class MapCounter<T> implements Serializable{

		private static final long serialVersionUID = 1L;

		Map<T, Integer> counts;

		public MapCounter() {
			counts = new HashMap<T, Integer>();
		}

		public MapCounter(Map<T, Integer> map) {
			counts = map;
		}


		public void increment(T key){
			if(!counts.containsKey(key))
				counts.put(key, 0);
			counts.put(key, counts.get(key)+1);
		}

		public void add(T key, Integer amount){
			if(!counts.containsKey(key))
				counts.put(key, 0);
			counts.put(key, counts.get(key)+amount);
		}

		public Map<T,Integer> getMap(){
			return counts;
		}

		public int get(T toCount){
			return counts.containsKey(toCount) ? counts.get(toCount) : 0;
		}

		public void set(T toCount,int i){
			counts.put(toCount,i);
		}

		public Set<T> keySet() {
			return counts.keySet();
		}

		public Collection<Integer> values() {
			return counts.values();
		}

		public void addAll(MapCounter<T> count) {
			for(T o : count.keySet()){
				add(o,count.get(o));
			}
		}

		public String toString() {
			return counts.toString();
		}

		public void clear() {
			counts.clear();
		}

		public void remove(Object o) {
			counts.remove(o);
		}
	}

	public static class MapSummer<T>{
		Map<T, Double> sum;

		public MapSummer() {
			sum = new HashMap<T, Double>();
		}

		public MapSummer(Map<T,Double> map) {
			sum = map;
		}

		public void add(T key, Double amount){
			if(!sum.containsKey(key))
				sum.put(key, 0.0);
			sum.put(key, sum.get(key)+amount);
		}

		public Map<T,Double> getMap(){
			return sum;
		}

		public double get(T toCount){
			return sum.containsKey(toCount) ? sum.get(toCount) : 0.0;
		}

		public void put(T toCount,double d){
			sum.put(toCount,d);
		}

		public Collection<T> keySet() {
			return sum.keySet();
		}

		public Collection<Double> values() {
			return sum.values();
		}

		public void addAll(MapSummer<T> summer) {
			for(T o : summer.keySet()){
				add(o,summer.get(o));
			}
		}

		public String toString() {
			return sum.toString();
		}

		public void clear() {
			sum.clear();
		}

		public void remove(T o) {
			sum.remove(o);
		}
	}

	public static class ComparableComparator<T extends Comparable<T>> implements Comparator<T>{
		//		<T extends Comparable<T>>
		public int compare(T o1, T o2) {
			return  o1.compareTo(o2);
		}

	}

	public static class NewtonsMethod{

		int defaultIterations=20;
		double defaultError=1E-5;
		Function f;
		Function fPrime;

		public NewtonsMethod(Function f, Function fPrime) {
			this.f=f;
			this.fPrime=fPrime;
		}

		public NewtonsMethod(Function f, double dx) {
			this.f=f;
			class DerivativeFunction implements Function{
				Function f;
				double dx;
				public DerivativeFunction(Function f, double dx) {
					this.f=f;
					this.dx=dx;
				}

				public double derivative(double x, double dx){
					return (f.eval(x)-f.eval(dx))/dx;
				}

				public double derivative(double x){
					return derivative(x, dx);
				}

				public double eval(double x) {
					return derivative(x);
				}
			}
			this.fPrime=new DerivativeFunction(f, dx);
		}

		public double root(){
			return root(0, defaultIterations, defaultError);
		}

		public double root(double init, int nIterations, double error){
			double root = init;
			int i=0;
			while(i<nIterations && Math.abs(f.eval(root))>error){
				root = root - (f.eval(root)/fPrime.eval(root));
				i++;
				//				System.out.println(root);
				//				System.out.println(f.eval(root));
			}
			return root;		
		}

		public interface Function{
			public double eval(double x);
		}
	}

	/**
	 * 
	 * @author sol
	 *
	 * @param <T>
	 */
	public static class MapExtreme<T>{
		Comparator<T> comparator;
		Map<Object, T> extremeValues;
		boolean isMin;

		public MapExtreme(Comparator<T> comparator, boolean isMin) {
			extremeValues = new HashMap<Object, T>();
			this.comparator=comparator;
			this.isMin=isMin;
		}

		public void put(Object key, T value){
			if(!extremeValues.containsKey(key))
				extremeValues.put(key, value);
			else{
				int comparison = comparator.compare(extremeValues.get(key), value);
				if(comparison < 0 && !isMin || comparison > 0 && isMin)
					extremeValues.put(key, value );
			}
		}

		public Map<Object,T> getExtremes(){
			return extremeValues;
		}

		public T getExtreme(Object obj){
			return extremeValues.get(obj);
		}
	}

	public static class ExtremeTracker<T>{
		Comparator<T> comparator;
		T maxValue;
		T minValue;

		public boolean hasExtrema;

		public ExtremeTracker(Comparator<T> comparator) {
			maxValue=null;
			this.comparator=comparator;
			hasExtrema=false;
		}

		@SuppressWarnings({ "rawtypes", "unchecked" })
		public ExtremeTracker() {
			maxValue=null;
			this.comparator=new ComparableComparator();
			hasExtrema=false;
		}

		public void put(T value){
			if(!hasExtrema){
				maxValue=value;
				minValue=value;
				hasExtrema = true;
			}
			else{
				{
					int comparison = comparator.compare(maxValue, value);
					if(comparison < 0 ){
						maxValue=value;
					}
				}
				{
					int comparison = comparator.compare(minValue, value);
					if(comparison > 0 ){
						minValue=value;
					}
				}
			}
		}

		public T getMax(){
			return maxValue;
		}

		public T getMin(){
			return minValue;
		}
	}

	/**
	 * 
	 * @author sol
	 *
	 * @param <S> The Object type
	 * @param <T> The Value type
	 */
	public static class ExtremeObjectTracker<S,T>{
		Comparator<T> comparator;
		//		S extreme;
		S maxExtreme;
		Set<S> maxObjects;
		S minExtreme;
		Set<S> minObjects;
		//		T extremeValue;
		T maxValue;
		T minValue;
		boolean hasExtrema;

		public ExtremeObjectTracker(Comparator<T> comparator) {
			//			extremeValue=null;
			this.comparator=comparator;
			maxObjects=new HashSet<S>();
			minObjects=new HashSet<S>();
			hasExtrema=false;
			//			this.isMin=isMin;
		}

		//		public void put(S key, T value){
		//			if(extremeValue==null){
		//				extreme=key;
		//				extremeValue=value;
		//			}
		//			else{
		//				int comparison = comparator.compare(extremeValue, value);
		//				if(comparison < 0 && !isMin || comparison > 0 && isMin){
		//					extreme=key;
		//					extremeValue=value;
		//				}
		//			}
		//		}

		public boolean hasExtrema(){
			return hasExtrema;
		}
		
		public void put(S key, T value){
			if(maxExtreme==null){
				hasExtrema = true;
				//				extreme=key;
				maxExtreme=key;
				minExtreme=key;
				maxObjects.add(key);
				minObjects.add(key);
				//				extremeValue=value;
				maxValue=value;
				minValue=value;
			}
			else{
				int max_comparison = comparator.compare(maxValue, value);
				if(max_comparison < 0){ 
					maxObjects.clear();
					maxObjects.add(key);
					maxExtreme=key;
					maxValue=value;
				}
				else if (max_comparison == 0){
					maxObjects.add(key);
				}
				int min_comparison = comparator.compare(minValue, value);
				if(min_comparison > 0){
					minObjects.clear();
					minObjects.add(key);
					minExtreme=key;
					minValue=value;
				}
				else if(min_comparison == 0){
					minObjects.add(key);
				}
			}
		}

		public S getMaxObject(){
			return maxExtreme;
		}

		public S getMinObject(){
			return minExtreme;
		}

		public Set<S> getMaxObjects(){
			return maxObjects;
		}

		public Set<S> getMinObjects(){
			return minObjects;
		}

		public T getMax(){
			return maxValue;
		}

		public T getMin(){
			return minValue;
		}
	}

	public static class MapList<S,T> implements Serializable{
		Map<S,List<T>> mapList;
		public MapList() {
			mapList = new HashMap<S, List<T>>();
		}

		public MapList(Map<S,List<T>> map) {
			mapList = map;
		}

		public void put(S key, T value){
			if(!mapList.containsKey(key))
				mapList.put(key, new LinkedList<T>());
			mapList.get(key).add(value);
		}

		public void putAll(S key, Collection<T> values){
			for(T value : values)
				put(key,value);
		}

		public List<T> get(S key){
			if(!mapList.containsKey(key))
				return new ArrayList<T>();
			return mapList.get(key);
		}

		public Map<S,List<T>> getMap(){
			return mapList;
		}

		public Set<S> keySet(){
			return mapList.keySet();
		}

		public Collection<List<T>> values(){
			return mapList.values();
		}

		public void addAll(Map<S,List<T>> toAdd){
			for(S s : toAdd.keySet())
				putAll(s,toAdd.get(s));
		}

		public void addAll(MapList<S,T> toAdd){
			addAll(toAdd.mapList);
		}

		public String toString(){
			return mapList.toString();
		}

		public void remove(S key) {
			mapList.remove(key);
		}

		public int size() {
			return mapList.size();
		}
	}

	public static class MapSet<S,T>{
		public Map<S,Set<T>> mapList;
		public MapSet() {
			mapList = new HashMap<S, Set<T>>();
		}

		public void put(S key, T value){
			if(!mapList.containsKey(key))
				mapList.put(key, new HashSet<T>());
			mapList.get(key).add(value);
		}

		public void putAll(S key, Collection<T> values){
			for(T value : values)
				put(key,value);
		}

		public Set<T> get(S key){
			if(!mapList.containsKey(key))
				return new HashSet<T>();
			return mapList.get(key);
		}

		public Map<S,Set<T>> getMap(){
			return mapList;
		}

		public Set<S> keySet(){
			return mapList.keySet();
		}

		public Collection<Set<T>> values(){
			return mapList.values();
		}

		public void addAll(Map<S,Set<T>> toAdd){
			for(S s : toAdd.keySet())
				putAll(s,toAdd.get(s));
		}

		public void addAll(MapSet<S,T> toAdd){
			addAll(toAdd.mapList);
		}

		public String toString(){
			return mapList.toString();
		}

		public void remove(String key) {
			mapList.remove(key);
		}

		public void remove(String key, String value) {
			mapList.get(key).remove(value);
		}

		public void removeAll(String key, Collection<T> values) {
			mapList.get(key).removeAll(values);
		}

		public int size() {
			return mapList.size();
		}
	}


	public static class MapFactory <S,T>{
		Map<S,T> map;
		Instantiator<T> i;

		public MapFactory(Instantiator<T> i) {
			map = new HashMap<S,T>();
			this.i = i;			
		}

		public T get(S s){
			if(!map.containsKey(s)){
				map.put(s, i.instantiate());
			}

			return map.get(s);
		}

		/**
		 * 
		 * @param s
		 * @param objects - optional arguments to be supplied to the instantiator
		 * @return
		 */
		public T get(S s, Object...objects){
			if(!map.containsKey(s)){
				map.put(s, i.instantiate(objects));
			}

			return map.get(s);
		}

		public Map<S,T> getMap(){
			return map;
		}

	}

	public static interface Instantiator<T>{

		public T instantiate(Object...objects);

	}

	public static class DefaultInstantiator<T> implements Instantiator<T>{

		Constructor<T> constructor;
		//		public DefaultInstantiator(Class<T> c) throws SecurityException, NoSuchMethodException {
		//			constructor = c.getConstructor();
		//		}

		public DefaultInstantiator(Class<T> c, Class<?>...types) throws SecurityException, NoSuchMethodException {
			constructor = c.getConstructor(types);
		}
		
		public T instantiate(Object...objects) {		
			try {
				return constructor.newInstance(objects);
			} catch (IllegalArgumentException e) {
				e.printStackTrace();
			} catch (InstantiationException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			} catch (InvocationTargetException e) {
				e.printStackTrace();
			}
			throw new IllegalArgumentException();
		}

	}
	/**
	 * 
	 * @param toFilter
	 * @param criteria
	 * @return a list of all items that satisfy the criteria
	 */
	public static <T> List<T> filter(Collection<T> toFilter, Criteria<T> criteria){
		List<T> filtered = new ArrayList<T>();
		for(T t : toFilter){
			if(criteria.satisfies(t))
				filtered.add(t);
		}
		return filtered;
	}

	public static interface Function<S,T> {
		public T evaluate(S s);
	}

	public static <S,T> List<T> evaluate(Collection<S> domain, Function<S,T> function){
		List<T> values = new ArrayList<T>(domain.size());

		for(S s : domain){
			values.add(function.evaluate(s));
		}

		return values;
	}

	public static <S,T> List<T> convolve(List<S> values, Function<List<S>,T> function, int bandWidth){
		List<T> convolution = new ArrayList<T>(values.size());

		for(int i=0; i<values.size()-bandWidth+1; i++){
			convolution.add(function.evaluate(values.subList(i, i+bandWidth)));
		}

		return convolution;
	}
	
	public static <S extends Number,T> List<Double> convolution(List<S> values, Function<List<S>,T> function, int bandWidth, double[] p){
		
		List<Double> convolution = new ArrayList<Double>(values.size());
		IntegralTable it = new IntegralTable(p);
		
		for(int i=0; i<values.size(); i++){
			double c = 0;
			for(int j=Math.max(i-bandWidth, 0); j<Math.min(values.size(), i+bandWidth+1); j++){
				System.out.printf("%d %d %d\n", i,j,j-i+bandWidth);
				c+= values.get(j).doubleValue()*p[j-i+bandWidth];
			}
			System.out.printf("%d,%d %d\n", Math.max(-(i-bandWidth), 0), Math.min(values.size()-i+bandWidth,p.length-1), p.length);
			System.out.println(it.sum(Math.max(-(i-bandWidth), 0), Math.min(values.size()-i+bandWidth,p.length-1)));
			c/=it.sum(Math.max(-(i-bandWidth), 0), Math.min(values.size()-i+bandWidth,p.length-1));
//			c=Util.median(values.subList(Math.max(i-bandWidth, 0), Math.min(values.size(), i+bandWidth+1)));
				convolution.add(c);
//			convolution.add(function.evaluate(values.subList(Math.max(i-bandWidth, 0), Math.min(values.size(), i+bandWidth+1))));
		}

		return convolution;
	}

	/**
	 * 
	 * @param toCheck
	 * @param criteria
	 * @return whether any item in the list satisfies the criteria
	 */
	public static <T> boolean anySatisfy(Iterable<T> toCheck, Criteria<T> criteria){
		for(T t : toCheck){
			if(criteria.satisfies(t))
				return true;
		}
		return false;
	}

	/**
	 * 
	 * @param toCheck
	 * @param criteria
	 * @return whether any item in the list satisfies the criteria
	 */
	public static <T> boolean allSatisfy(Collection<T> toCheck, Criteria<T> criteria){
		for(T t : toCheck){
			if(!criteria.satisfies(t))
				return false;
		}
		return true;
	}

	public static interface Criteria<T>{
		public boolean satisfies(T t);
	}

	public static <S extends Number> List<Integer> convertNumber(Collection<S> numbers, int i) {
		List<Integer> converted = new ArrayList<Integer>(numbers.size());
		for(Number n : numbers){	
			converted.add(n.intValue());
		}
		return converted;
	}

	public static <S extends Number> List<Byte> convertNumber(Collection<S> numbers, byte b) {
		List<Byte> converted = new ArrayList<Byte>(numbers.size());
		for(Number n : numbers){	
			converted.add(n.byteValue());
		}
		return converted;
	}

	public static <S extends Number> List<Double> convertNumber(Collection<S> numbers, double b) {
		List<Double> converted = new ArrayList<Double>(numbers.size());
		for(Number n : numbers){	
			converted.add(n.doubleValue());
		}
		return converted;
	}

	public static <S extends Number> List<Short> convertNumber(Collection<S> numbers, short s) {
		List<Short> converted = new ArrayList<Short>(numbers.size());
		for(Number n : numbers){	
			converted.add(n.shortValue());
		}
		return converted;
	}

	public static <S extends Number> List<Long> convertNumber(Collection<S> numbers, long l) {
		List<Long> converted = new ArrayList<Long>(numbers.size());
		for(Number n : numbers){	
			converted.add(n.longValue());
		}
		return converted;
	}

	public static <S extends Number> List<Float> convertNumber(Collection<S> numbers, float f) {
		List<Float> converted = new ArrayList<Float>(numbers.size());
		for(Number n : numbers){	
			converted.add(n.floatValue());
		}
		return converted;
	}

	public static List<Color> colorGradient(int n) {
		return colorGradient(n,1);
	}

	public static List<Color> colorGradient(int n, double a) {
		List<Color> gradient = new ArrayList<Color>(n);
		for(int i=0; i < n; i++){
			gradient.add(hsba(i*1./n, 1, 1, a));
		}
		return gradient;
	}

	public static String wrap(String string, int width) {
		int nNewLines = (string.length()+width-1)/width;
		StringBuffer wrappedString = new StringBuffer(string.length()+nNewLines);
		for (int i = 0; i < nNewLines; i++) {
			wrappedString.append(string.subSequence(i*width, Math.min((i*width)+width,string.length())));
			wrappedString.append("\n");
		}
		wrappedString.deleteCharAt(wrappedString.length()-1);
		return wrappedString.toString();
	}

	public static void copy(File src, File dest) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(src));
		BufferedWriter bw = new BufferedWriter(new FileWriter(dest));
		int c = br.read();
		while(c!=-1){
			bw.write(c);
			c=br.read();
		}
		br.close();
		bw.close();
	}

	public static <S,T> List<T> bin(List<S> list, int binWidth, Function<List<S>, T> binningFunction) {
		int l=(list.size()+binWidth-1)/binWidth;

		List<T> binnedList = new ArrayList<T>(l);
		for(int i=0; i<l; i++){
			binnedList.add(binningFunction.evaluate(list.subList(i*binWidth, Math.min((i*binWidth) + binWidth, list.size()))));
		}

		return binnedList;
	}

	public static <T> void initialize(T[] array, Instantiator<T> instantiator, Object... objects) {
		for(int i=0; i<array.length; i++){
			array[i] = instantiator.instantiate(objects);
		}
	}


}
