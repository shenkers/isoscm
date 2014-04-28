package tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

import util.IO;
import util.Util;
import util.Util.Function;

//chr1    4836814 4836815 7.206|AATAAA    118     +
public class ParseBed {
	public static class BEDIterator implements Iterable<AnnotatedRegion>, Iterator<AnnotatedRegion>{
		Scanner scan;
		IO.LineTokenizer tokenizer;
		boolean hasNext;

		public BEDIterator(File bedFile, boolean ignoreHeader) throws FileNotFoundException {
			scan = IO.bufferedScanner(bedFile);
			tokenizer = new IO.LineTokenizer(scan, "\t");
			if(ignoreHeader)
				tokenizer.next();
			hasNext=true;
		}

		public BEDIterator(File bedFile) throws FileNotFoundException {
			scan = IO.bufferedScanner(bedFile);
			tokenizer = new IO.LineTokenizer(scan, "\t");
			hasNext=true;
		}

		public BEDIterator(String bedFile) throws FileNotFoundException {
			this(new File(bedFile));
		}

		public boolean hasNext() {
			return hasNext && tokenizer.hasNext();
		}

		public AnnotatedRegion next() {
			String[] item = tokenizer.next();

			String chr = item[0];
			String annotation = item.length > 3 ? item[3] : null;
			int start = Integer.parseInt(item[1])+1;
			int end = Integer.parseInt(item[2]);
			char strand = item.length > 5 ? item[5].charAt(0) : '.';

			AnnotatedRegion region = new AnnotatedRegion(annotation, chr, start, end, strand);

			if(!tokenizer.hasNext()){
				hasNext=false;
				scan.close();
			}

			return region;
		}

		public void remove() {

		}

		public Iterator<AnnotatedRegion> iterator() {
			return this;
		}

	}
	
//	Iterable<List<AnnotatedRegion>>,
	public static class BED12Iterable implements Iterable<List<AnnotatedRegion>>{
		String bedFile;
		
		public BED12Iterable(String bedFile) {
			this.bedFile = bedFile;
		}
		
		public Iterator<List<AnnotatedRegion>> iterator() {
			try {
				return new BED12Iterator(bedFile);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			return null;
		}
		
	}
	
	public static class BED12Iterator implements Iterator<List<AnnotatedRegion>>{
		Scanner scan;
		IO.LineTokenizer tokenizer;
		boolean hasNext;

		Function<String, Integer> parseInt;

		public BED12Iterator(File bedFile) throws FileNotFoundException {
			scan = IO.bufferedScanner(bedFile);
			tokenizer = new IO.LineTokenizer(scan, "\t");
			parseInt = new Function<String, Integer>() {
				public Integer evaluate(String s) {
					return Integer.parseInt(s);
				}
			};
			hasNext=true;
		}

		public BED12Iterator(String bedFile) throws FileNotFoundException {
			this(new File(bedFile));
		}

		public boolean hasNext() {
			return hasNext && tokenizer.hasNext();
		}

		public List<AnnotatedRegion> next() {
			String[] item = tokenizer.next();

			String chr = item[0];
			String annotation = item.length > 3 ? item[3] : null;
			int start = Integer.parseInt(item[1])+1;
			int end = Integer.parseInt(item[2]);
			char strand = item.length > 5 ? item[5].charAt(0) : '.';

			List<Integer> lengths = Util.evaluate(Util.list(item[10].split(",")), parseInt);
			List<Integer> offsets = Util.evaluate(Util.list(item[11].split(",")), parseInt);
						
			List<AnnotatedRegion> exons = new LinkedList<AnnotatedRegion>();
			for(int i=0; i<lengths.size(); i++){
				exons.add(new AnnotatedRegion(annotation, chr, start+offsets.get(i), start+offsets.get(i)+lengths.get(i)-1, strand));
			}

			if(!tokenizer.hasNext()){
				hasNext=false;
				scan.close();
			}

			return exons;
		}

		public void remove() {

		}
	}

	public static void main(String[] args) throws FileNotFoundException {
		BEDIterator gi = new BEDIterator("/home/sol/data/derti/processed_files/mouse/GSM747481_mouse_brain.sites.clustered.bed");
		while(gi.hasNext)
			System.out.println(gi.next());
	}
}
