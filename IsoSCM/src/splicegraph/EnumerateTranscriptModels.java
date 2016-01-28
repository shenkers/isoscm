package splicegraph;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import processing.FindSpliceJunctions;
import tools.AnnotatedRegion;
import tools.BEDTools.BEDWriter;
import tools.GTFTools.GTFWriter;
import tools.GenomicIntervalSet;
import tools.ParseBed.BEDIterator;
import tools.ParseGTF.TranscriptIterator;
import tools.IntervalTools;
import tools.StrandedGenomicIntervalSet;
import tools.StrandedGenomicIntervalTree;
import util.Util;
import util.Util.MapCounter;

public class EnumerateTranscriptModels {


	public static void filterRegion(String chr, int start, int end, String segmentationBed) throws FileNotFoundException{
		BEDIterator bi = new BEDIterator("test3.bed");
		BEDWriter bw = new BEDWriter(segmentationBed);
		for(AnnotatedRegion r : bi){
			if(r.chr.equals(chr) && r.start >= start && r.end <= end){
				System.out.println(r);
				bw.write(r.chr, r.start, r.end, '.');
			}
		}
		bw.close();
	}

	public static void filterJunctions(SamReader sfr, String chr, int start, int end, String junctionBed) throws FileNotFoundException{
		List<AnnotatedRegion> junctions = FindSpliceJunctions.spliceJunction(sfr, chr, start, end, true);

		BEDWriter bw = new BEDWriter(junctionBed);
		for(AnnotatedRegion r : junctions){
			bw.write(r.chr, r.start, r.end, r.strand);
		}
		bw.close();
	}

	public static void writeTranscript(StrandedGenomicIntervalSet exonIntervals, BEDWriter bw){
		String chr = null;
		Character strand = null;
		int start = Integer.MAX_VALUE;
		int end = Integer.MIN_VALUE;

		List<Integer> starts = new LinkedList<Integer>();
		List<Integer> lengths = new LinkedList<Integer>();
		for(AnnotatedRegion r : exonIntervals){
			start = Math.min(start,r.start);
			end = Math.max(end,r.end);
			starts.add(r.start);
			lengths.add(r.size);
			strand = r.strand;
			chr = r.chr;
		}

		bw.write("", chr, start, end, starts, lengths, strand);
	}

	//	public static void writeTranscript(StrandedGenomicIntervalSet exonIntervals, List<Integer> segmentBoundaries, BEDWriter bw){
	//		for(Collection<Integer> boundary : new Util.subsetIterator<Integer>(segmentBoundaries, 2)){
	//		int min = Collections.min(boundary);
	//		int max = Collections.max(boundary);
	//		
	//		StrandedGenomicIntervalSet exonIntervals = new Strand
	//		
	//		writeTranscript(exonIntervals, bw);		
	//	}
	//	}

	public static double fraction_using_jnct(SamReader sfr, String chr, int start, int end){
		SAMRecordIterator sri = sfr.query(chr, start, start, false);
		Set<String> splice_junction_reads = new HashSet<String>();
		Set<String> non_splice_junction_reads = new HashSet<String>();

		int n_splice = 0; 
		int n_non = 0;

		while(sri.hasNext()){
			SAMRecord sr = sri.next();

			if(sr.getAlignmentStart() == start)
				continue;

			String readID = sr.getAlignmentStart()+sr.getCigarString();

			Cigar cigar = sr.getCigar();
			int alignmentPosition = sr.getAlignmentStart();

			boolean usesSpliceJnct = false;

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
					if(isIntron){
						if(alignmentPosition==start && alignmentPosition+cigarElement.getLength()-1 == end){
							usesSpliceJnct = true;
						}
					}
					alignmentPosition += cigarElement.getLength();
				}			
			}

			if(usesSpliceJnct){
				n_splice++;
				//				splice_junction_reads.add(readID);	
			}
			else{
				n_non++;
				//				non_splice_junction_reads.add(readID);
			}
		}
		sri.close();

		sri = sfr.query(chr, end, end, false);

		while(sri.hasNext()){
			SAMRecord sr = sri.next();

			if(sr.getAlignmentEnd() == end)
				continue;

			String readID = sr.getAlignmentStart()+sr.getCigarString();
			Cigar cigar = sr.getCigar();
			int alignmentPosition = sr.getAlignmentStart();

			boolean usesSpliceJnct = false;

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
					if(isIntron){
						if(alignmentPosition==start && alignmentPosition+cigarElement.getLength()-1 == end){
							usesSpliceJnct = true;
						}
					}
					alignmentPosition += cigarElement.getLength();
				}			
			}

			if(usesSpliceJnct){
				//				splice_junction_reads.add(readID);	
			}
			else{
				n_non++;
				//				non_splice_junction_reads.add(readID);
			}
		}
		sri.close();

		//		return splice_junction_reads.size()*1./(splice_junction_reads.size()+non_splice_junction_reads.size());
		return n_splice*1./(n_splice+n_non);
	}

	public static void enumerateSpliceModels(String chr, int start, int end, List<AnnotatedRegion> splices, BEDWriter bw) throws FileNotFoundException{

		//		BEDWriter bw = new BEDWriter("/dev/stdout");
		for(char strand : Util.list('+','-')){
			List<AnnotatedRegion> spliceJunctions = new LinkedList<AnnotatedRegion>();
			for(AnnotatedRegion j : splices){
				if(j.strand==strand && j.start > start && j.end < end)
					spliceJunctions.add(j);
			}

			Collections.sort(spliceJunctions, new Comparator<AnnotatedRegion>() {

				public int compare(AnnotatedRegion arg0, AnnotatedRegion arg1) {
					if(arg0.start!=arg1.start)
						return arg0.start-arg1.start;
					else
						return arg0.end-arg1.end;
				}
			});

			StrandedGenomicIntervalSet exonModel = new StrandedGenomicIntervalSet();
			exonModel.add(chr, start, end, strand);

			//			writeTranscript(exonModel, bw);
			if(spliceJunctions.size()>0)
				enumerateSpliceModels(chr, start, end, strand, 0, spliceJunctions, exonModel, bw);

			//		bw.write(name, chr, start, end, exon_starts, exon_lengths, strand)
		}

	}

	public static void enumerateSpliceModels(String chr, int start, int end, char strand, List<AnnotatedRegion> splices, BEDWriter bw) throws FileNotFoundException{

		//		BEDWriter bw = new BEDWriter("/dev/stdout");
		List<AnnotatedRegion> spliceJunctions = new LinkedList<AnnotatedRegion>();
		for(AnnotatedRegion j : splices){
			if(j.strand==strand && j.start > start && j.end < end)
				spliceJunctions.add(j);
		}

		Collections.sort(spliceJunctions, new Comparator<AnnotatedRegion>() {

			public int compare(AnnotatedRegion arg0, AnnotatedRegion arg1) {
				if(arg0.start!=arg1.start)
					return arg0.start-arg1.start;
				else
					return arg0.end-arg1.end;
			}
		});

		StrandedGenomicIntervalSet exonModel = new StrandedGenomicIntervalSet();
		exonModel.add(chr, start, end, strand);

		//			writeTranscript(exonModel, bw);
		if(spliceJunctions.size()>0)
			enumerateSpliceModels(chr, start, end, strand, 0, spliceJunctions, exonModel, bw);

		//		bw.write(name, chr, start, end, exon_starts, exon_lengths, strand)

	}

	public static void enumerateSpliceModels(String chr, int start, int end, char strand, int spliceIndex, List<AnnotatedRegion> splices, StrandedGenomicIntervalSet exonModel, BEDWriter bw) throws FileNotFoundException{

		if(spliceIndex == splices.size()){
			writeTranscript(exonModel, bw);		
		} 
		else{
			AnnotatedRegion junction = null;
			int i=spliceIndex;
			for(; i<splices.size(); i++){
				AnnotatedRegion j = splices.get(i);
				if(j.start > start && j.end < end){
					junction = j;
					break;
				}
			}

			if(junction==null){
				writeTranscript(exonModel, bw);	
			}
			else{
				exonModel.remove(junction.chr, junction.start, junction.end, junction.strand);
				enumerateSpliceModels(chr, junction.end+1, end, strand, i+1, splices, exonModel, bw);


				boolean conflicts = false;

				// conflicts if its end is between another start and end
				if(conflicts)
					exonModel.add(junction.chr, junction.start, junction.end, junction.strand);

				int conflictEnd = junction.end;

				for(int j=i+1; j<splices.size() && splices.get(j).start < conflictEnd && splices.get(j).end < end; j++){
					//					System.out.println("CONFLICT");
					exonModel.add(junction.chr, junction.start, junction.end, junction.strand);
					junction = splices.get(j);
					conflictEnd = Math.min(conflictEnd,junction.end);
					exonModel.remove(junction.chr, junction.start, junction.end, junction.strand);
					enumerateSpliceModels(chr, junction.end+1, end, strand, j+1, splices, exonModel, bw);
				}

				//				System.out.printf("esm\t%s\n", junction);
			}
		}
	}

	public static void old_main(String[] args) throws FileNotFoundException {
		//		chr1:192,805,797-192,858,788

		//		intervals("test3.gtf");
		SamReader sfr = SamReaderFactory.makeDefault().open(new File("/home/sol/data/sangercenter/hippocampus.bam"));

		//		filterRegion("chr1", 192805797, 192858788);
		//		filterJunctions(sfr,"chr1", 192805797, 192858788);

		BEDIterator junction_it = new BEDIterator("assembly_junctions.bed");

		StrandedGenomicIntervalTree<Map<String,Object>> junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();
		List<AnnotatedRegion> splices = new LinkedList<AnnotatedRegion>();
		MapCounter<String> mc = new MapCounter<String>();
		for(AnnotatedRegion junction : junction_it){
			boolean contained = false;
			for(AnnotatedRegion overlap : junctions.overlappingRegions(junction.chr, junction.start, junction.end, junction.strand)){
				if(overlap.start == junction.start && overlap.end == junction.end){
					contained=true;
				}
			}
			if(!contained){
				System.out.println(junction);
				junctions.add(junction);
				splices.add(junction);
			}
			mc.increment(junction.toString());
		}

		Collections.sort(splices, new Comparator<AnnotatedRegion>() {

			public int compare(AnnotatedRegion arg0, AnnotatedRegion arg1) {
				if(arg0.start!=arg1.start)
					return arg0.start-arg1.start;
				else
					return arg0.end-arg1.end;
			}
		});

		Iterator<AnnotatedRegion> si = splices.iterator();
		while(si.hasNext()){
			AnnotatedRegion next = si.next();
			String nextS = next.toString();
			int n = mc.get(nextS);
			double usage = fraction_using_jnct(sfr, next.chr, next.start, next.end);

			if(n<3 || usage < .15)
				si.remove();


			System.out.printf("%s\t%d\t%.2f\n", nextS,mc.get(nextS),usage);
		}
		//		System.exit(0);


		List<Integer> segmentBoundaries = new LinkedList<Integer>();
		GenomicIntervalSet gis = new GenomicIntervalSet();
		BEDIterator bi = new BEDIterator("assembly_segmentation.bed");
		for(AnnotatedRegion r : bi){
			gis.add(r.chr, r.start, r.start);
			gis.add(r.chr, r.end, r.end);
		}

		for(AnnotatedRegion r : splices){
			gis.remove(r.chr, r.start-50, r.start+50);
			gis.remove(r.chr, r.end-50, r.end+50);
		}

		System.out.println("boundaries");
		for(AnnotatedRegion r : gis){
			segmentBoundaries.add(r.start);
		}
		Collections.sort(segmentBoundaries);

		BEDWriter bw = new BEDWriter("assembly.bed");
		for(Collection<Integer> boundary : new Util.subsetIterator<Integer>(segmentBoundaries, 2)){
			int min = Collections.min(boundary);
			int max = Collections.max(boundary);

			enumerateSpliceModels("chr1", min, max, splices, bw);
		}

		//		enumerateSpliceModels("chr1", 192805797, 192858788, splices, segmentBoundaries, bw);

		bw.close();
		//		for(AnnotatedRegion r : splices){
		//			System.out.println(r);
		//		}
	}


	public static void main(String[] args) throws FileNotFoundException {

	
		String loc = "chr1:183,269,781-183,413,660";
		loc="chr1:134,869,918-134,928,273";
		loc="chr1:120,176,951-120,213,058";
		loc="chr1:95,650,693-95,689,757";
		loc="chr1:95,685,317-95,699,849";
		loc="chr1:88,235,609-88,268,664";
		loc="chr1:80,249,083-80,347,157";
		loc = loc.replaceAll(",", "");

		String chr = loc.split(":")[0];
		int start = Integer.parseInt(loc.split(":")[1].split("-")[0]);
		int end = Integer.parseInt(loc.split(":")[1].split("-")[1]);

		String segmentationBed = "test_segmentation.bed";
		String junctionBed = "test_junctions.bed";
		//		chr1:192,805,797-192,858,788

		//		intervals("test3.gtf");
		SamReader sfr = SamReaderFactory.makeDefault().open(new File("/home/sol/data/sangercenter/hippocampus.bam"));

		filterRegion(chr, start, end, segmentationBed);
		filterJunctions(sfr,chr, start, end, junctionBed);

		BEDIterator junction_it = new BEDIterator(junctionBed);

		StrandedGenomicIntervalTree<Map<String,Object>> junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();
		List<AnnotatedRegion> splices = new LinkedList<AnnotatedRegion>();
		MapCounter<String> mc = new MapCounter<String>();
		for(AnnotatedRegion junction : junction_it){
			boolean contained = false;
			for(AnnotatedRegion overlap : junctions.overlappingRegions(junction.chr, junction.start, junction.end, junction.strand)){
				if(overlap.start == junction.start && overlap.end == junction.end){
					contained=true;
				}
			}
			if(!contained){
				System.out.println(junction);
				junctions.add(junction);
				splices.add(junction);
			}
			mc.increment(junction.toString());
		}

		Collections.sort(splices, new Comparator<AnnotatedRegion>() {

			public int compare(AnnotatedRegion arg0, AnnotatedRegion arg1) {
				if(arg0.start!=arg1.start)
					return arg0.start-arg1.start;
				else
					return arg0.end-arg1.end;
			}
		});

		Iterator<AnnotatedRegion> si = splices.iterator();
		while(si.hasNext()){
			AnnotatedRegion next = si.next();
			String nextS = next.toString();
			int n = mc.get(nextS);
			double usage = fraction_using_jnct(sfr, next.chr, next.start, next.end);

			if(n<3 || (n<20&&usage < .15))
				si.remove();


			System.out.printf("%s\t%d\t%.2f\n", nextS,mc.get(nextS),usage);
		}
		//		System.exit(0);


		List<Integer> segmentBoundaries = new LinkedList<Integer>();
		GenomicIntervalSet gis = new GenomicIntervalSet();
		BEDIterator bi = new BEDIterator(segmentationBed);
		for(AnnotatedRegion r : bi){
			gis.add(r.chr, r.start, r.start);
			gis.add(r.chr, r.end, r.end);
		}

		for(AnnotatedRegion r : splices){
			gis.remove(r.chr, r.start-50, r.start+50);
			gis.remove(r.chr, r.end-50, r.end+50);
		}

		for(AnnotatedRegion r : gis){
			segmentBoundaries.add(r.start);
		}
		Collections.sort(segmentBoundaries);

		String output = "test_assembly.bed";

		BEDWriter bw = new BEDWriter(output);
		for(Collection<Integer> boundary : new Util.subsetIterator<Integer>(segmentBoundaries, 2)){
			int min = Collections.min(boundary);
			int max = Collections.max(boundary);

			enumerateSpliceModels(chr, min, max, splices, bw);
		}

		//		enumerateSpliceModels("chr1", 192805797, 192858788, splices, segmentBoundaries, bw);

		bw.close();
		//		for(AnnotatedRegion r : splices){
		//			System.out.println(r);
		//		}
	}
}
