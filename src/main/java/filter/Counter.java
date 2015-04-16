package filter;

public class Counter<T> implements Consumer<T> {

	public int n;

	public Counter() {
		n=0;
	}

	public void add(T t) {
		n++;
	}
	
	public void reset(){
		n=0;
	}

}
