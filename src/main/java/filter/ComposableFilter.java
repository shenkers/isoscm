package filter;

import java.util.LinkedList;
import java.util.List;

public abstract class ComposableFilter<T> implements Consumer<T>, Emitter<T>{

	List<Consumer<T>> listeners;
	
	public ComposableFilter() {
		listeners = new LinkedList<Consumer<T>>();
	}
	
	/**
	 * if t satisfies this filter, process it, then
	 * send it to all the listeners
	 * @param t
	 */
	public void add(T t){
		if(satisfies(t)){
			emit(t);
		}
	}
	
	/**
	 * 
	 * @param t
	 * @return true if t passes the filter
	 */
	public abstract boolean satisfies (T t);
	
	public void emit(T t){
		for(Consumer<T> cf : listeners){
			cf.add(t);
		}
	}
	
	public void append(Consumer<T> c){
		listeners.add(c);
	}
}
