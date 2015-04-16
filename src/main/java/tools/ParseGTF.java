package tools;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.IO;
import util.Util;


public class ParseGTF {

	public static void parseGTF(Collection<String> roi) throws FileNotFoundException{
		PrintStream out = IO.bufferedPrintstream("all_gene_intron_properties.txt");
		Scanner scan = IO.bufferedScanner("/Users/sol/Desktop/apa_tracks/ginormous.100806.gtf");

		Map<String,Map<String,List<AnnotatedRegion>>> geneTranscriptModels = new HashMap<String, Map<String,List<AnnotatedRegion>>>();

		for(String[] item : new IO.LineTokenizer(scan, "\t")){
			//System.out.println(item[2]);
			//System.out.println(Util.list(item));
			//	System.out.println(Util.list(item[0],item[2],item[3],item[4]));
			//	System.out.println(Util.list(item[8].split(";")));

			String geneName = null;
			String transcriptName = null;

			for(String it : Util.list(item[8].split(";\\W*"))){
				int startI = it.indexOf("\"")+1;
				int endI = it.lastIndexOf("\"");
				if(it.matches("^gene_id.*")){
					geneName= it.substring(startI,endI);
				}
				if(it.matches("^transcript_id.*")){
					transcriptName = it.substring(startI,endI);
				}
			}

			if(!roi.contains(geneName))
				continue;

			char strand = item[6].equals("+") ? '+' : '-';

			if(!geneTranscriptModels.containsKey(geneName))
				geneTranscriptModels.put(geneName,new HashMap<String, List<AnnotatedRegion>>());
			if(!geneTranscriptModels.get(geneName).containsKey(transcriptName))
				geneTranscriptModels.get(geneName).put(transcriptName, new ArrayList<AnnotatedRegion>());

			AnnotatedRegion ar = new AnnotatedRegion(item[2], item[0], Integer.parseInt(item[3]), Integer.parseInt(item[4]), strand);
			geneTranscriptModels.get(geneName).get(transcriptName).add(ar);
			//		if(transcriptName.equals("CG31687.c"))
			//			System.out.println(ar);
		}
		//System.out.println(geneTranscriptModels);
		scan.close();
		out.println("transcript\tUTRSize\tn-Introns\tminSize\tmaxSize\tmeanSize");
		for(String gene : geneTranscriptModels.keySet()){
			//for(String transcript : geneTranscriptModels.get(gene).keySet()){
			String transcript = longestUtrIsoform(geneTranscriptModels.get(gene));
			AnnotatedRegion UTR3 = getAnnotatedRegion(geneTranscriptModels.get(gene).get(transcript), "3UTR");
			List<AnnotatedRegion> introns = getIntrons(geneTranscriptModels.get(gene).get(transcript));
			int nIntrons = introns.size();
			List<Integer> sizes = new ArrayList<Integer>();
			for(AnnotatedRegion intron : introns){
				sizes.add(intron.size);
			}
			int minSize = sizes.size() > 0 ? Collections.min(sizes) : 0;
			int maxSize = sizes.size() > 0 ? Collections.max(sizes) : 0;
			double meanSize = sizes.size() > 0 ? Util.sum(sizes)/sizes.size() : 0;
			out.printf("%s\t%d\t%d\t%d\t%d\t%.1f\n",transcript,UTR3==null ? 0 : UTR3.size,nIntrons,minSize,maxSize,meanSize);
			//	}
		}
		out.close();
	}

	public static Map<String,Map<String,List<AnnotatedRegion>>> getTranscripts(String gtfFile, Collection<String> roi) throws FileNotFoundException{
		//		PrintStream out = IO.bufferedPrintstream("all_gene_intron_properties.txt");
		Scanner scan = IO.bufferedScanner(gtfFile);

		Map<String,Map<String,List<AnnotatedRegion>>> geneTranscriptModels = new HashMap<String, Map<String,List<AnnotatedRegion>>>();

		for(String[] item : new IO.LineTokenizer(scan, "\t")){
			//System.out.println(item[2]);
			//System.out.println(Util.list(item));
			//	System.out.println(Util.list(item[0],item[2],item[3],item[4]));
			//	System.out.println(Util.list(item[8].split(";")));

			String geneName = null;
			String transcriptName = null;

			for(String it : Util.list(item[8].split(";\\W*"))){
				int startI = it.indexOf("\"")+1;
				int endI = it.lastIndexOf("\"");
				if(it.matches("^gene_id.*")){
					geneName= it.substring(startI,endI);
				}
				if(it.matches("^transcript_id.*")){
					transcriptName = it.substring(startI,endI);
				}
			}

			if(!roi.contains(geneName) && !roi.contains(transcriptName))
				continue;

			char strand = item[6].equals("+") ? '+' : '-';

			if(!geneTranscriptModels.containsKey(geneName))
				geneTranscriptModels.put(geneName,new HashMap<String, List<AnnotatedRegion>>());
			if(!geneTranscriptModels.get(geneName).containsKey(transcriptName))
				geneTranscriptModels.get(geneName).put(transcriptName, new ArrayList<AnnotatedRegion>());

			AnnotatedRegion ar = new AnnotatedRegion(item[2], item[0], Integer.parseInt(item[3]), Integer.parseInt(item[4]), strand);
			geneTranscriptModels.get(geneName).get(transcriptName).add(ar);
			//		if(transcriptName.equals("CG31687.c"))
			//			System.out.println(ar);
		}
		//System.out.println(geneTranscriptModels);
		scan.close();
		//		out.println("transcript\tUTRSize\tn-Introns\tminSize\tmaxSize\tmeanSize");
		//		for(String gene : geneTranscriptModels.keySet()){
		//			//for(String transcript : geneTranscriptModels.get(gene).keySet()){
		//			String transcript = longestUtrIsoform(geneTranscriptModels.get(gene));
		//			AnnotatedRegion UTR3 = getAnnotatedRegion(geneTranscriptModels.get(gene).get(transcript), "3UTR");
		//			List<AnnotatedRegion> introns = getIntrons(geneTranscriptModels.get(gene).get(transcript));
		//			int nIntrons = introns.size();
		//			List<Integer> sizes = new ArrayList<Integer>();
		//			for(AnnotatedRegion intron : introns){
		//				sizes.add(intron.size);
		//			}
		//			int minSize = sizes.size() > 0 ? Collections.min(sizes) : 0;
		//			int maxSize = sizes.size() > 0 ? Collections.max(sizes) : 0;
		//			double meanSize = sizes.size() > 0 ? Util.sum(sizes)/sizes.size() : 0;
		////			out.printf("%s\t%d\t%d\t%d\t%d\t%.1f\n",transcript,UTR3==null ? 0 : UTR3.size,nIntrons,minSize,maxSize,meanSize);
		//			//	}
		//		}
		//		out.close();

		return geneTranscriptModels;
	}

	public static Map<String,Map<String,List<AnnotatedRegion>>> getTranscripts(String gtfFile) throws FileNotFoundException{
		Scanner scan = IO.bufferedScanner(gtfFile);

		Map<String,Map<String,List<AnnotatedRegion>>> geneTranscriptModels = new HashMap<String, Map<String,List<AnnotatedRegion>>>();

		for(String[] item : new IO.LineTokenizer(scan, "\t")){
			String geneName = null;
			String transcriptName = null;

			for(String it : Util.list(item[8].split(";\\W*"))){
				int startI = it.indexOf("\"")+1;
				int endI = it.lastIndexOf("\"");
				if(it.matches("^gene_id.*")){
					geneName= it.substring(startI,endI);
				}
				if(it.matches("^transcript_id.*")){
					transcriptName = it.substring(startI,endI);
				}
			}

			char strand = item[6].equals("+") ? '+' : '-';

			if(!geneTranscriptModels.containsKey(geneName))
				geneTranscriptModels.put(geneName,new HashMap<String, List<AnnotatedRegion>>());
			if(!geneTranscriptModels.get(geneName).containsKey(transcriptName))
				geneTranscriptModels.get(geneName).put(transcriptName, new ArrayList<AnnotatedRegion>());

			AnnotatedRegion ar = new AnnotatedRegion(item[2], item[0], Integer.parseInt(item[3]), Integer.parseInt(item[4]), strand);
			geneTranscriptModels.get(geneName).get(transcriptName).add(ar);
		}

		scan.close();

		return geneTranscriptModels;
	}

	public static Map<String,Map<String,List<AnnotatedRegion>>> getFlybaseTranscripts(String gtfFile, Collection<String> roi) throws FileNotFoundException{
		//		PrintStream out = IO.bufferedPrintstream("all_gene_intron_properties.txt");
		Scanner scan = IO.bufferedScanner(gtfFile);

		Map<String,Map<String,List<AnnotatedRegion>>> geneTranscriptModels = new HashMap<String, Map<String,List<AnnotatedRegion>>>();

		for(String[] item : new IO.LineTokenizer(scan, "\t")){
			//System.out.println(item[2]);
			//System.out.println(Util.list(item));
			//	System.out.println(Util.list(item[0],item[2],item[3],item[4]));
			//	System.out.println(Util.list(item[8].split(";")));

			String geneName = null;
			String transcriptName = null;

			for(String it : Util.list(item[8].split(";\\W*"))){
				int startI = it.indexOf("\"")+1;
				int endI = it.lastIndexOf("\"");
				if(it.matches("^gene_id.*")){
					geneName= it.substring(startI,it.lastIndexOf('-'));
				}
				if(it.matches("^transcript_id.*")){
					transcriptName = it.substring(startI,endI);
				}
			}

			if(!roi.contains(geneName) && !roi.contains(transcriptName))
				continue;

			char strand = item[6].equals("+") ? '+' : '-';

			if(!geneTranscriptModels.containsKey(geneName))
				geneTranscriptModels.put(geneName,new HashMap<String, List<AnnotatedRegion>>());
			if(!geneTranscriptModels.get(geneName).containsKey(transcriptName))
				geneTranscriptModels.get(geneName).put(transcriptName, new ArrayList<AnnotatedRegion>());

			AnnotatedRegion ar = new AnnotatedRegion(item[2], item[0], Integer.parseInt(item[3]), Integer.parseInt(item[4]), strand);
			geneTranscriptModels.get(geneName).get(transcriptName).add(ar);
			//		if(transcriptName.equals("CG31687.c"))
			//			System.out.println(ar);
		}
		//System.out.println(geneTranscriptModels);
		scan.close();
		//		out.println("transcript\tUTRSize\tn-Introns\tminSize\tmaxSize\tmeanSize");
		//		for(String gene : geneTranscriptModels.keySet()){
		//			//for(String transcript : geneTranscriptModels.get(gene).keySet()){
		//			String transcript = longestUtrIsoform(geneTranscriptModels.get(gene));
		//			AnnotatedRegion UTR3 = getAnnotatedRegion(geneTranscriptModels.get(gene).get(transcript), "3UTR");
		//			List<AnnotatedRegion> introns = getIntrons(geneTranscriptModels.get(gene).get(transcript));
		//			int nIntrons = introns.size();
		//			List<Integer> sizes = new ArrayList<Integer>();
		//			for(AnnotatedRegion intron : introns){
		//				sizes.add(intron.size);
		//			}
		//			int minSize = sizes.size() > 0 ? Collections.min(sizes) : 0;
		//			int maxSize = sizes.size() > 0 ? Collections.max(sizes) : 0;
		//			double meanSize = sizes.size() > 0 ? Util.sum(sizes)/sizes.size() : 0;
		////			out.printf("%s\t%d\t%d\t%d\t%d\t%.1f\n",transcript,UTR3==null ? 0 : UTR3.size,nIntrons,minSize,maxSize,meanSize);
		//			//	}
		//		}
		//		out.close();

		return geneTranscriptModels;
	}

	public static class TranscriptIterator implements Iterator<AnnotatedRegion>, Iterable<AnnotatedRegion> {
		Scanner scan;
		IO.LineTokenizer tokenizer;
		private boolean hasNext;
		Pattern attributePattern1 = Pattern.compile("([^ ]+) \"([^\"]*)\"");
		Pattern attributePattern2 = Pattern.compile("([^=]+)=(.+)");

		public TranscriptIterator(String gtfFile) throws FileNotFoundException{
			scan = IO.bufferedScanner(gtfFile);
			tokenizer = new IO.LineTokenizer(scan, "\t");
			hasNext = true;
		}
		
		public TranscriptIterator(File gtfFile) throws FileNotFoundException{
			scan = IO.bufferedScanner(gtfFile);
			tokenizer = new IO.LineTokenizer(scan, "\t");
			hasNext = true;
		}


		public boolean hasNext() {
			return hasNext && tokenizer.hasNext();
		}

		public AnnotatedRegion next() {
			String[] item = tokenizer.next();

			Map<String,Object> attributes = new HashMap<String, Object>();
			for(String it : Util.list(item[8].split(";\\W*"))){
				Matcher m1 = attributePattern1.matcher(it);
				Matcher m2 = attributePattern2.matcher(it);
				if(m1.find()){
					attributes.put(m1.group(1), m1.group(2));
				}
				else if(m2.find()){
					attributes.put(m2.group(1), m2.group(2));
				}

			}

			char strand = item[6].charAt(0);


			//	System.out.println(Util.list(item));
			AnnotatedRegion ar = new AnnotatedRegion(item[2], item[0], Integer.parseInt(item[3]), Integer.parseInt(item[4]), strand, attributes);

			if(!tokenizer.hasNext()){
				hasNext = false;
				scan.close();	
			}
			return ar;
		}

		public void remove() {
			// DO NOTHING
		}


		public Iterator<AnnotatedRegion> iterator() {
			return this;
		}
	}

	public static void parseOtherGTF() throws FileNotFoundException{
		PrintStream out = IO.bufferedPrintstream("all_gene_intron_properties.txt");
		Scanner scan = IO.bufferedScanner("/Users/sol/Desktop/apa_tracks/ginormous.100806.gtf");

		out.println("transcript\tUTRSize\tn-Introns\tminSize\tmaxSize\tmeanSize");

		String lastTranscriptName=null;
		List<AnnotatedRegion> model = new ArrayList<AnnotatedRegion>();
		for(String[] item : new IO.LineTokenizer(scan, "\t")){
			//System.out.println(item[2]);
			//System.out.println(Util.list(item));
			//	System.out.println(Util.list(item[0],item[2],item[3],item[4]));
			//	System.out.println(Util.list(item[8].split(";")));

			String geneName = null;
			String transcriptName = null;

			for(String it : Util.list(item[8].split(";\\W*"))){
				int startI = it.indexOf("\"")+1;
				int endI = it.lastIndexOf("\"");
				if(it.matches("^gene_id.*")){
					geneName= it.substring(startI,endI);
				}
				if(it.matches("^transcript_id.*")){
					transcriptName = it.substring(startI,endI);
				}
			}

			if(lastTranscriptName==null)
				lastTranscriptName=transcriptName;

			char strand = item[6].equals("+") ? '+' : '-';

			AnnotatedRegion ar = new AnnotatedRegion(item[2], item[0], Integer.parseInt(item[3]), Integer.parseInt(item[4]), strand);

			if(!transcriptName.equals(lastTranscriptName)){
				AnnotatedRegion UTR3 = getAnnotatedRegion(model, "3UTR");
				List<AnnotatedRegion> introns = getIntrons(model);
				int nIntrons = introns.size();
				List<Integer> sizes = new ArrayList<Integer>();
				for(AnnotatedRegion intron : introns){
					sizes.add(intron.size);
				}
				int minSize = sizes.size() > 0 ? Collections.min(sizes) : 0;
				int maxSize = sizes.size() > 0 ? Collections.max(sizes) : 0;
				double meanSize = sizes.size() > 0 ? Util.sum(sizes)/sizes.size() : 0;
				out.printf("%s\t%d\t%d\t%d\t%d\t%.1f\n",transcriptName,UTR3==null ? 0 : UTR3.size,nIntrons,minSize,maxSize,meanSize);
				model.clear();

			}

			lastTranscriptName = transcriptName;
			model.add(ar);
			//		if(transcriptName.equals("CG31687.c"))
			//			System.out.println(ar);
		}
		//System.out.println(geneTranscriptModels);
		scan.close();
		out.close();
	}

	public static AnnotatedRegion getAnnotatedRegion(List<AnnotatedRegion> model, String annotation){
		AnnotatedRegion region = null;
		for(AnnotatedRegion r  : model){
			if(r.annotation.equals(annotation))
				region=r;
		}
		return region;
	}

	private static String longestUtrIsoform(Map<String, List<AnnotatedRegion>> isoforms) {
		List<String> names = Util.list(isoforms.keySet());
		List<List<AnnotatedRegion>> annotations = Util.get(isoforms, names);
		String longest = names.get(annotations.indexOf(Collections.max(annotations, new Comparator<List<AnnotatedRegion>>() {

			public int compare(List<AnnotatedRegion> iso1, List<AnnotatedRegion> iso2) {
				AnnotatedRegion UTR1 = null; AnnotatedRegion UTR2 = null;
				UTR1 = getAnnotatedRegion(iso1, "3UTR");
				UTR2 = getAnnotatedRegion(iso2, "3UTR");
				if(UTR1==null && UTR2==null)
					return 0;
				else if(UTR1==null ^ UTR2==null)
					return (UTR1==null ? -1 : 1);
				else
					return UTR1.size-UTR2.size;
			}

		})));
		return longest;
	}

	public static int[] getStartEndIndexes(List<AnnotatedRegion> model){
		int minIndex = Integer.MAX_VALUE;
		int maxIndex = Integer.MIN_VALUE;
		for(AnnotatedRegion ar : model){
			minIndex = Collections.min(Util.list(minIndex,ar.start,ar.end));
			maxIndex = Collections.max(Util.list(minIndex,ar.start,ar.end));
		}
		return new int[]{minIndex,maxIndex};
	}

	public static List<AnnotatedRegion> getIntrons(List<AnnotatedRegion> model){
		int[] startEndIndexes = getStartEndIndexes(model);
		char strand = model.iterator().next().strand;
		List<AnnotatedRegion> introns = new ArrayList<AnnotatedRegion>();

		List<Integer> indexes = new ArrayList<Integer>();
		for(AnnotatedRegion r : model){
			if(r.annotation.equals("5UTR") || r.annotation.equals("3UTR") || r.annotation.equals("CDS") || r.annotation.equals("stop_codon")){
				indexes.add(r.start); indexes.add(r.end);
			}
		}
		Collections.sort(indexes);
		for(int index=2; index < indexes.size(); index = index + 2){
			int intronSize = indexes.get(index) - indexes.get(index-1) - 1;
			if(intronSize > 0){
				introns.add(new AnnotatedRegion("intron",model.get(0).chr,indexes.get(index-1)+1,indexes.get(index)-1, strand));
			}
		}

		return introns;
	}

	public static List<String> genes() throws FileNotFoundException{
		List<String> genes = new ArrayList<String>();
		Scanner scan = IO.bufferedScanner("genes_of_interest");
		while(scan.hasNextLine()){
			String line = scan.nextLine();
			genes.add(line);
		}
		return genes;
	}

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException {
		List<String> genes = genes();
		System.out.println(genes);
		//parseGTF(genes);
		parseOtherGTF();

	}

	public boolean hasNext() {
		// TODO Auto-generated method stub
		return false;
	}

	public Map<String, Map<String, List<AnnotatedRegion>>> next() {
		// TODO Auto-generated method stub
		return null;
	}

	public void remove() {
		// TODO Auto-generated method stub

	}

}


