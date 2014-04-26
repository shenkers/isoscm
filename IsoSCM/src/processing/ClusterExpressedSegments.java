package processing;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.NormalDistribution;

import cern.colt.list.DoubleArrayList;
import cern.jet.random.engine.RandomEngine;
import cern.jet.stat.quantile.DoubleQuantileFinder;
import cern.jet.stat.quantile.QuantileFinderFactory;
import filter.ComposableFilter;
import filter.ConsumingReadFilter;
import filter.Counter;
import filter.MultiMappingfilter;
import filter.SAMRecordStrandednessFilter;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import splicegraph.ExonSpliceGraph.ConnectedComponentResult;
import tools.AnnotatedRegion;
import tools.BAMTools;
import tools.BEDTools.BEDWriter;
import tools.GTFTools.GTFWriter;
import tools.GenomicIntervalSet;
import tools.GenomicIntervalTree;
import tools.GenomicPositionTree;
import tools.IntervalTools;
import tools.StrandedSegmentClusterer;
import tools.ParseBed.BEDIterator;
import tools.ParseGTF.TranscriptIterator;
import tools.PositionClustering.StrandedClusterer;
import tools.StrandedGenomicIntervalSet;
import tools.StrandedGenomicIntervalTree;
import tools.StrandedGenomicPositionTree;
import tools.Strandedness;
import util.ChartUtils;
import util.IO;
import util.Util;
import util.Stats.MomentTracker;
import util.Util.ExtremeObjectTracker;
import util.Util.ExtremeTracker;
import util.Util.MapList;

public class ClusterExpressedSegments {

	public static int nConsumingReads(SAMFileReader sfr, Strandedness strandedness, String chr, int start, int end, boolean contained, boolean isNegativeStrand){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			SAMRecord sr = sri.next();

			Boolean isNegStrand = null;

			switch (strandedness) {
			case reverse_forward:
				isNegStrand = sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag() ? false : true;
				break;

			default:
				throw new RuntimeException("Unimplimented");
			}

			if(!(isNegativeStrand ^ isNegStrand)){

				int alignmentPosition = sr.getAlignmentStart();
				Cigar cigar = sr.getCigar();

				boolean consumesOverlappingBases = false;


				consumesOverlappingBases:{
					for(int i=0; i<cigar.numCigarElements(); i++){
						CigarElement cigarElement = cigar.getCigarElement(i);
						if(cigarElement.getOperator().consumesReferenceBases()){
							boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
							for(int j=0; j<cigarElement.getLength(); j++){
								if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
									consumesOverlappingBases=true;
									break consumesOverlappingBases;
								}
								alignmentPosition++;
							}
						}			
					}
				}

				if(consumesOverlappingBases)
					n++;

			}
		}

		sri.close();
		return n;
	}

	public static GenomicIntervalTree<Map<String,Object>> intervals(String gtf, GTFWriter gw) throws FileNotFoundException {
		StrandedGenomicIntervalTree<Map<String,Object>> git = new StrandedGenomicIntervalTree<Map<String,Object>>();
		GenomicIntervalTree<Map<String,Object>> git2 = new GenomicIntervalTree<Map<String,Object>>();

		TranscriptIterator ti = new TranscriptIterator(gtf);

		for(AnnotatedRegion r : ti){
			if(!git.contains(r.chr, r.start, r.end,r.strand)){
				git.add(r.chr, r.start, r.end,r.strand, r.attributes);
			}
		}

		// remove all stranded changepoints that share one boundary 
		// with an unstranded segment
		// and
		// remove all stranded changepoints that are shorter than
		// another segment sharing the same boundary
		List<AnnotatedRegion> irrelevant = new LinkedList<AnnotatedRegion>();
		for(AnnotatedRegion r : git){
			if(r.strand=='.'){
				for(char c : Util.list('+','-')){
					int r5p = r.get5Prime(c);
					for(AnnotatedRegion n : git.overlappingRegions(r.chr, r5p, r5p, c)){
						if(r5p == n.get5Prime()){
							irrelevant.add(n);
						}
					}
				}
			}
			else{
				int r5p = r.get5Prime();
				for(AnnotatedRegion n : git.overlappingRegions(r.chr, r5p, r5p, r.strand)){
					if(r5p == n.get5Prime() && r.size < n.size){
						irrelevant.add(r);
					}
				}
			}
		}

		for(AnnotatedRegion r : irrelevant){
			git.remove(r.chr, r.start, r.end, r.strand);
		}

		// find those that are on both strands,
		// delete them, those that remain are those that are only on one strand.
		// for each of these, assign them to the closest up/down-stream segment boundary.

		StrandedGenomicPositionTree unmatchedSegmentBoundaries = new StrandedGenomicPositionTree();

		for(AnnotatedRegion r : git){
			if(r.strand=='.'){
				unmatchedSegmentBoundaries.add(r.chr, r.start, '+');
				unmatchedSegmentBoundaries.add(r.chr, r.end, '-');
			}
			else{
				unmatchedSegmentBoundaries.add(r.chr, r.get5Prime(), r.strand);
			}
		}

		GenomicIntervalTree<Map<String,Object>> matchedSegmentBoundaries = new GenomicIntervalTree<Map<String,Object>>();
		for(AnnotatedRegion r : git){
			if(r.strand!='.'){
				AnnotatedRegion closest = unmatchedSegmentBoundaries.getClosestUpstream(r.chr, r.get3Prime(), r.isNegativeStrand()?'+':'-');
				//			System.out.println(r);
				//			System.out.printf("\t%s\n", closest);
				if(closest!=null){

					String chr = r.chr;
					int start = Math.min(r.start,closest.start);
					int end = Math.max(r.end,closest.end);

					boolean nested = false;
					boolean matched = false;
					for(AnnotatedRegion matchedSegment : git.overlappingRegions(chr, start, end, r.isNegativeStrand()?'+':'-')){
						if(matchedSegment.get5Prime()==closest.start){
							matched |= true;
							//							System.out.println(matchedSegment);
							//							System.out.println(matchedSegment.toAttributeString());
							closest.addAttribute("expression", matchedSegment.getAttribute("expression"));
							int inter_5p_overlap = Math.abs(r.get5Prime() - matchedSegment.get5Prime());
							if(inter_5p_overlap < matchedSegment.size){
								System.out.printf("nested %s in %s\n", r, matchedSegment);
								nested = true;
							}

						}
					}
					if(matched && !(nested  || matchedSegmentBoundaries.contains(chr, start, end))){
						//						System.out.printf("(%s,%s)\t%s:%d-%d\t%s\t%s\n",r,closest,r.chr,Math.min(r.start,closest.start), Math.max(r.end,closest.end),r.toAttributeString(),closest.toAttributeString());
						Map<String, Object> attributes = new HashMap<String, Object>();
						attributes.put("consolidated", true);
						if(r.isNegativeStrand()){
							attributes.put("start_expression", closest.getAttribute("expression"));
							attributes.put("end_expression", r.getAttribute("expression"));
						}
						else{
							attributes.put("start_expression", r.getAttribute("expression"));
							attributes.put("end_expression", closest.getAttribute("expression"));
						}
						matchedSegmentBoundaries.add(r.chr, Math.min(r.start,closest.start), Math.max(r.end,closest.end), attributes);
					}
				}
			}
		}
		//		System.exit(0);

		for(AnnotatedRegion r : matchedSegmentBoundaries){
			//			System.out.println(r);
			gw.write("exon", r.chr, r.start, r.end, '.', r.toGTFAttributeString());
		}

		//		BEDWriter bw = new BEDWriter("test3.bed");
		//		for(AnnotatedRegion r : git2){
		//			bw.write(r.chr, r.start, r.end, '.');
		//		}
		//		bw.close();


		return git2;
	}

	static class ConsolidatedIntervalResult{
		GenomicIntervalTree<Map<String,Object>> segments;
		GenomicIntervalTree<Map<String,Object>> changePoints;
		GenomicIntervalTree<Map<String,Object>> increasePoints;
		GenomicIntervalTree<Map<String,Object>> decreasePoints;

		public ConsolidatedIntervalResult(GenomicIntervalTree<Map<String, Object>> segments, GenomicIntervalTree<Map<String, Object>> changePoints, GenomicIntervalTree<Map<String, Object>> increasePoints, GenomicIntervalTree<Map<String, Object>> decreasePoints) {
			this.segments = segments;
			this.changePoints = changePoints;
			this.increasePoints = increasePoints;
			this.decreasePoints = decreasePoints;
		}
	}

	public static ConsolidatedIntervalResult expressedSegmentIntervals(String gtf) throws FileNotFoundException {
		StrandedGenomicIntervalTree<Map<String,Object>> git = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> unique_git = new StrandedGenomicIntervalTree<Map<String,Object>>();
		GenomicPositionTree gpt = new GenomicPositionTree();
		GenomicIntervalTree<Map<String,Object>> max_git = new GenomicIntervalTree<Map<String,Object>>();
		GenomicIntervalTree<Map<String,Object>> increase_git = new GenomicIntervalTree<Map<String,Object>>();
		GenomicIntervalTree<Map<String,Object>> decrease_git = new GenomicIntervalTree<Map<String,Object>>();
		GenomicIntervalTree<Map<String,Object>> expressedSegmentIntervals = new GenomicIntervalTree<Map<String,Object>>();

		TranscriptIterator ti = new TranscriptIterator(gtf);

		for(AnnotatedRegion r : ti){
			if(!unique_git.contains(r.chr, r.start, r.end,r.strand)){
				unique_git.add(r.chr, r.start, r.end,r.strand);
				gpt.add(r.chr, r.start);
			}
			git.add(r.chr, r.start, r.end,r.strand, r.attributes);
		}

		for(AnnotatedRegion r : unique_git){
			MapList<String, Double> mles = new MapList<String, Double>();

			//			System.out.println(r);
			for(AnnotatedRegion g : git.overlappingRegions(r.chr, r.start, r.end, r.strand)){
				for(String s : Util.list("before_mle","after_mle")){				
					mles.put(s, Double.parseDouble((String) g.getAttribute(s)));					
				}
			}

			Map<String,Object> attributes = new HashMap<String, Object>();
			for(String s : Util.list("before_mle","after_mle")){
				attributes.put(s, Collections.max(mles.get(s)));
			}	

			boolean isIncrease = (Double) attributes.get("after_mle") > (Double) attributes.get("before_mle");


			if(!max_git.contains(r.chr, r.start, r.start)){
				max_git.add(r.chr,r.start,r.end,attributes);
				if(isIncrease)
					increase_git.add(r.chr,r.start,r.end,attributes);
				else
					decrease_git.add(r.chr,r.start,r.end,attributes);
				//				gw.write("exon", r.chr, r.start, r.end, '.', AnnotatedRegion.GTFAttributeString(attributes));		
			}
		}

		//		System.out.println();

		for(AnnotatedRegion r : max_git){
			//			System.out.println(r);

			AnnotatedRegion upstream = gpt.getClosestUpstream(r.chr, r.start-1);
			//			AnnotatedRegion downstream = gpt.getClosestDownstream(r.chr, r.end+1);

			if(upstream!=null){
				Map<String,Object> attributes = new HashMap<String, Object>();
				attributes.put("mle", Math.max((Double) r.getAttribute("before_mle"), (Double) max_git.overlappingRegions(upstream.chr, upstream.start, upstream.end, upstream.strand).iterator().next().getAttribute("after_mle")));		
				//				gw.write("exon", r.chr, upstream.start, r.start, '.', AnnotatedRegion.GTFAttributeString(attributes));		
				expressedSegmentIntervals.add(r.chr, upstream.start, r.start, attributes);
			}

			//			if(downstream!=null){
			//				Map<String,Object> attributes = new HashMap<String, Object>();
			//				attributes.put("mle", Math.max((Double) r.getAttribute("after_mle"), (Double) max_git.overlappingRegions(downstream.chr, downstream.start, downstream.end, downstream.strand).iterator().next().getAttribute("before_mle")));
			//				gw.write("exon", r.chr, r.start, downstream.start, '.', AnnotatedRegion.GTFAttributeString(attributes));		
			//			
			//			}


			//			System.out.printf("%s\t%s\t%s\n", gpt.getClosestUpstream(r.chr, r.start-1),r,gpt.getClosestDownstream(r.chr, r.end+1));
		}

		ConsolidatedIntervalResult result = new ConsolidatedIntervalResult(expressedSegmentIntervals, max_git, increase_git, decrease_git);

		return result;
	}

	public static void enumerateExons(String chr, int start, int end, char strand, List<AnnotatedRegion> splices, BEDWriter bw) throws FileNotFoundException{

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
			enumerateExons(chr, start, end, strand, 0, spliceJunctions, exonModel, bw);

		//		bw.write(name, chr, start, end, exon_starts, exon_lengths, strand)

	}

	public static void enumerateExons(String chr, int start, int end, char strand, int spliceIndex, List<AnnotatedRegion> splices, StrandedGenomicIntervalSet exonModel, BEDWriter bw) throws FileNotFoundException{

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
				enumerateExons(chr, junction.end+1, end, strand, i+1, splices, exonModel, bw);


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
					enumerateExons(chr, junction.end+1, end, strand, j+1, splices, exonModel, bw);
				}

				//				System.out.printf("esm\t%s\n", junction);
			}
		}
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

	public static StrandedGenomicIntervalTree<Map<String,Object>> readBedFileAsIntervalTree(File bedFile) throws FileNotFoundException{
		StrandedGenomicIntervalTree<Map<String,Object>> intervals = new StrandedGenomicIntervalTree<Map<String,Object>>(); 

		BEDIterator bi = new BEDIterator(bedFile);
		for(AnnotatedRegion r : bi){
			intervals.add(r);
		}

		return intervals;
	}

	public static boolean intervalContains(int position, int start, int end){
		return position >= start && position <= end;
	}

	public static void inferSegmentStrand(File acc_jnct_gtf, File segment_bed, BEDWriter inferred_exons) throws FileNotFoundException{
		TranscriptIterator sjs = new TranscriptIterator(acc_jnct_gtf);
		BEDIterator segs = new BEDIterator(segment_bed);
		StrandedGenomicIntervalTree<Map<String,Object>> sj = IntervalTools.buildRegionsTree(sjs, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> seg = IntervalTools.buildRegionsTree(segs, true, false);

		StrandedGenomicIntervalTree<Map<String,Object>> inferred = new StrandedGenomicIntervalTree<Map<String,Object>>();
		for(AnnotatedRegion r : sj){

			for(AnnotatedRegion e : seg.overlappingRegions(r.chr, r.start-1, r.start-1, '.')){
				if(!inferred.contains(e.chr, e.start, e.end, r.strand)){
					inferred.add(e.chr,e.start,e.end,r.strand);
					//					System.out.printf("%s %c\n", e, r.strand);
					inferred_exons.write(e.chr, e.start, e.end, r.strand);
				}
			}
			for(AnnotatedRegion e : seg.overlappingRegions(r.chr, r.end+1, r.end+1, '.')){
				if(!inferred.contains(e.chr, e.start, e.end, r.strand)){
					inferred.add(e.chr,e.start,e.end,r.strand);
					//					System.out.printf("%s %c\n", e, r.strand);
					inferred_exons.write(e.chr, e.start, e.end, r.strand);
				}
			}
		}

		// add all segments not in inferred already as unstranded

		for(AnnotatedRegion e : seg){
			boolean inferred_strand = false;
			for(char strand : Util.list('+','-')){
				inferred_strand |= inferred.contains(e.chr, e.start, e.end, strand);
			}
			if(!inferred_strand){
				inferred.add(e.chr,e.start,e.end,e.strand);
				//				System.out.printf("%s %c\n", e, e.strand);
				inferred_exons.write(e.chr, e.start, e.end, e.strand);
			}
		}

	}

	public static void trimInferredExons(File spliced_exon_gtf, GTFWriter trimmed_exons) throws FileNotFoundException{
		TranscriptIterator ti = new TranscriptIterator(spliced_exon_gtf);
		StrandedGenomicIntervalTree<Map<String,Object>> e = IntervalTools.buildRegionsTree(ti, false, true);
		StrandedGenomicIntervalSet regions = IntervalTools.buildStrandedIntervalSet(e);

		StrandedGenomicIntervalTree<Map<String,Object>> t5p = IntervalTools.buildTerminiTree(e, true, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> t3p = IntervalTools.buildTerminiTree(e, false, true, false);

		GenomicIntervalSet antisense = new GenomicIntervalSet();
		for(AnnotatedRegion p : regions){
			if(p.strand=='+'){
				for(AnnotatedRegion m : regions.overlappingRegions(p.chr, p.start, p.end, '-')){		
					int[] i = IntervalTools.intersection(p.start, p.end, m.start, m.end);
					antisense.add(p.chr, i[0], i[1]);
				}
			}

		}
		//		System.exit(0);

		for(AnnotatedRegion x : e){

			AnnotatedRegion trimmed = null;
			// if the exon needs to be trimmed
			if(antisense.overlappingRegions(x.chr, x.start, x.end).iterator().hasNext()){

				// for each exon in the region on either strand
				// figure out what type it is
				//					System.out.println(x.toAttributeString());
				String type = (String) x.getAttribute("type");
				if("3p_exon".equals(type)){
					AnnotatedRegion c = t5p.getClosestUpstream(x.chr, x.get5Prime(), IntervalTools.opposite(x.strand));
					trimmed = IntervalTools.makeRegion(x.chr, x.get5Prime(), c.start, x.isNegativeStrand());
					//												System.out.printf("before %s\n", x);
					//											System.out.printf("trimmed %s\n", trimmed);
				}
				else if("5p_exon".equals(type)){
					AnnotatedRegion c = t3p.getClosestDownstream(x.chr, x.get3Prime(), IntervalTools.opposite(x.strand));
					trimmed = IntervalTools.makeRegion(x.chr, c.start, x.get3Prime(), x.isNegativeStrand());
					//					trimmed_exons.write("exon", trimmed.chr, trimmed.start, trimmed.end, trimmed.strand, AnnotatedRegion.GTFAttributeString(x.attributes));
					//							
				}
				else if("internal_exon".equals(type)){
					//					System.out.printf("before %s\n", x);
					//						System.out.printf("trimmed %s\n", trimmed);
					trimmed = x;
				}

			}
			else{
				// write the exon unmodified
				trimmed = x;
			}

			trimmed_exons.write("exon", trimmed.chr, trimmed.start, trimmed.end, trimmed.strand, AnnotatedRegion.GTFAttributeString(x.attributes));

		}
	}

	public static void identifySplicedExons(File acc_jnct_gtf, File segment_bed, GTFWriter exons) throws FileNotFoundException{

		StrandedGenomicIntervalTree<Map<String, Object>> expressedSegments = readBedFileAsIntervalTree(segment_bed);

		StrandedGenomicIntervalTree<Map<String,Object>> sj = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> splice5ps = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> splice3ps = new StrandedGenomicIntervalTree<Map<String,Object>>();

		TranscriptIterator splice_junctions = new TranscriptIterator(acc_jnct_gtf);
		for(AnnotatedRegion splice_junction : splice_junctions){
			if(!sj.contains(splice_junction.chr, splice_junction.start, splice_junction.end, splice_junction.strand))
				sj.add(splice_junction);
			if(!splice5ps.contains(splice_junction.chr, splice_junction.get5Prime(), splice_junction.get5Prime(), splice_junction.strand))
				splice5ps.add(splice_junction.chr, splice_junction.get5Prime(), splice_junction.get5Prime(), splice_junction.strand);
			if(!splice3ps.contains(splice_junction.chr, splice_junction.get3Prime(), splice_junction.get3Prime(), splice_junction.strand))
				splice3ps.add(splice_junction.chr, splice_junction.get3Prime(), splice_junction.get3Prime(), splice_junction.strand);
		}

		for(AnnotatedRegion r : expressedSegments){
			int i= r.isNegativeStrand()?-1:1;
			StrandedGenomicIntervalTree<Map<String,Object>> written = new StrandedGenomicIntervalTree<Map<String,Object>>();

			for(AnnotatedRegion splice5p : splice5ps.overlappingRegions(r.chr, r.start-1, r.end+1, r.strand)){
				AnnotatedRegion splice3p = splice3ps.getClosestUpstream(splice5p.chr, splice5p.start, r.strand);
				if(splice3p!=null && intervalContains(splice3p.start, r.start-1, r.end+1)){
					while(splice3p!=null && intervalContains(splice3p.start, r.start-1, r.end+1)){
						// write internal exon
						int start = Math.min(splice5p.start,splice3p.start);
						int end = Math.max(splice5p.start,splice3p.start); 
						if(end-start>1 && !written.contains(r.chr, start+1, end-1, r.strand)){
							Map<String,Object> attributes = new HashMap<String, Object>();
							attributes.put("type", "internal_exon");
							// as long as there is no nested splice junction
							if(IntervalTools.ContainedIntervals(sj, r.chr, start, end, r.strand).size()==0){
								written.add(r.chr, start+1, end-1, r.strand);
								exons.write("exon", r.chr, start+1, end-1, r.strand,AnnotatedRegion.GTFAttributeString(attributes));
							}
						}
						splice3p = splice3ps.getClosestUpstream(r.chr, splice3p.start-i, r.strand);
					}
				}
				else{
					// write 5' exon				
					int start = Math.min(splice5p.start+1,r.get5Prime());
					int end = Math.max(splice5p.start-1,r.get5Prime());
					if(end-start>-1 && !written.contains(r.chr, start, end, r.strand)){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("type", "5p_exon");
						// as long as there is no nested splice junction
						if(IntervalTools.ContainedIntervals(sj, r.chr, start, end, r.strand).size()==0){
							written.add(r.chr, start, end, r.strand);
							exons.write("exon",r.chr, start, end, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
						}
					}
				}
			}

			for(AnnotatedRegion splice3p : splice3ps.overlappingRegions(r.chr, r.start-1, r.end+1, r.strand)){
				AnnotatedRegion splice5p = splice5ps.getClosestDownstream(splice3p.chr, splice3p.start, r.strand);
				if(splice5p!=null && intervalContains(splice5p.start, r.start-1, r.end+1)){
					while(splice5p!=null && intervalContains(splice5p.start, r.start-1, r.end+1)){
						// write internal exon
						int start = Math.min(splice5p.start,splice3p.start);
						int end = Math.max(splice5p.start,splice3p.start); 
						if(end-start>1 && !written.contains(r.chr, start+1, end-1, r.strand)){
							Map<String,Object> attributes = new HashMap<String, Object>();
							attributes.put("type", "internal_exon");
							// as long as there is no nested splice junction
							if(IntervalTools.ContainedIntervals(sj, r.chr, start, end, r.strand).size()==0){
								written.add(r.chr, start+1, end-1, r.strand);
								exons.write("exon", r.chr, start+1, end-1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
						}
						splice5p = splice5ps.getClosestDownstream(r.chr, splice5p.start+i, r.strand);
					}
				}
				else{
					// write 3' exon
					int start = Math.min(splice3p.start+1,r.get3Prime());
					int end = Math.max(splice3p.start-1,r.get3Prime()); 
					if(end-start>-1 && !written.contains(r.chr, start, end, r.strand)){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("type", "3p_exon");
						// as long as there is no nested splice junction
						if(IntervalTools.ContainedIntervals(sj, r.chr, start, end, r.strand).size()==0){
							written.add(r.chr, start, end, r.strand);
							exons.write("exon",r.chr, start, end, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
						}
					}
				}
			}
		}
	}

	public static void mergeSegments(File segment_bed, File merge_regions, BEDWriter merged_segments, int merge_radius) throws FileNotFoundException {
		StrandedSegmentClusterer ssc = new StrandedSegmentClusterer(merge_radius);
		BEDIterator bi = new BEDIterator(segment_bed);

		// cluster all segments together that are less than the merge radius apart 
		for(AnnotatedRegion r : bi){
			ssc.addSegment(r.chr, r.start, r.end, r.strand);
		}

		if(merge_regions!=null){


			StrandedSegmentClusterer merge_clusters = new StrandedSegmentClusterer(merge_radius);
			BEDIterator mr = new BEDIterator(merge_regions);
			for(AnnotatedRegion r : mr){
				merge_clusters.addSegment(r.chr, r.start, r.end, '.');
			}

			// if regions to merge across are specified, identify gaps that can be explained by known regions
			for(AnnotatedRegion r : merge_clusters.clusters){
				for(char strand : Util.list('+','-','.')){
					ExtremeTracker<Integer> starts = new ExtremeTracker<Integer>();
					ExtremeTracker<Integer> ends = new ExtremeTracker<Integer>();
					for(AnnotatedRegion segment : ssc.clusters.overlappingRegions(r.chr, r.start, r.end, strand)){
						starts.put(segment.start);
						ends.put(segment.end);
					}

					if(starts.hasExtrema && starts.getMax() > ends.getMin() && ends.getMin() - starts.getMax() < 3*merge_radius ){
						ssc.addSegment(r.chr, ends.getMin(), starts.getMax(), strand, false);
					}

				}
			}
		}

		for(AnnotatedRegion r : ssc){
			merged_segments.write(r.chr, r.start, r.end, r.strand);
		}

		merged_segments.close();
	}

	public static void identifyIntronicExons(File segment_bed, File acc_jnct_gtf, File spliced_exon_gtf, SAMFileReader bam, Strandedness strandedness, GTFWriter intronic_exon_writer, int minExtensionLength, double minExtensionFraction) throws FileNotFoundException{

		StrandedGenomicIntervalTree<Map<String, Object>> expressedSegments = readBedFileAsIntervalTree(segment_bed);
		StrandedGenomicIntervalTree<Map<String, Object>> sj = IntervalTools.buildRegionsTree(new TranscriptIterator(acc_jnct_gtf), true, false);
		StrandedGenomicIntervalTree<Map<String, Object>> splicedExons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> splicedExonBoundaries = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> written = new StrandedGenomicIntervalTree<Map<String,Object>>();

		TranscriptIterator ti = new TranscriptIterator(spliced_exon_gtf.getAbsolutePath());

		for(AnnotatedRegion spliced_exon : ti){
			if("internal_exon".equals(spliced_exon.getAttribute("type"))){
				splicedExons.add(spliced_exon);
				splicedExonBoundaries.add(spliced_exon.chr, spliced_exon.start, spliced_exon.start, spliced_exon.strand);
				splicedExonBoundaries.add(spliced_exon.chr, spliced_exon.end, spliced_exon.end, spliced_exon.strand);
			}
		}



		for(AnnotatedRegion segment : expressedSegments){
			AnnotatedRegion closestStart = splicedExonBoundaries.getClosestDownstream(segment.chr, segment.get5Prime(), segment.strand);
			AnnotatedRegion closestEnd = splicedExonBoundaries.getClosestUpstream(segment.chr, segment.get3Prime(), segment.strand);

			int i = segment.isNegativeStrand() ? 1 : -1;

			if(closestStart != null && Math.abs(closestStart.start - segment.get5Prime()) > minExtensionLength && intervalContains(closestStart.start, segment.start, segment.end) && !splicedExons.overlappingRegions(segment.chr, segment.get5Prime(), segment.get5Prime(), segment.strand).iterator().hasNext()){

				// the extension read density is taken to be the density from the closest exon border to the boundary of 
				// the extension, with an extra border given by minExtensionLength
				int s = Math.min(segment.get5Prime(), closestStart.start+(i*minExtensionLength)); 
				int e = Math.max(segment.get5Prime(), closestStart.start+(i*minExtensionLength));
				//	a			double extension_read_density = nConsumingReads(bam, strandedness, segment.chr, s, e, false, segment.isNegativeStrand())*1./(e-s+1);
				for(AnnotatedRegion spliced_exon : splicedExons.overlappingRegions(segment.chr, segment.start, segment.end, segment.strand)){
					int start = Math.min(segment.get5Prime(), spliced_exon.get3Prime());
					int end = Math.max(segment.get5Prime(), spliced_exon.get3Prime());

					if(IntervalTools.ContainedIntervals(sj, segment.chr, start, end, segment.strand).size()==0 && !written.contains(segment.chr, start, end, segment.strand)){
						//	a					double exon_read_density = nConsumingReads(bam, strandedness, spliced_exon.chr, spliced_exon.start, spliced_exon.end, false, segment.isNegativeStrand())*1./spliced_exon.size;
						//	a					if(extension_read_density/exon_read_density > minExtensionFraction){
						//							System.out.printf("%s:%d-%d\n", spliced_exon.chr, Math.min(segment.get5Prime(), spliced_exon.get5Prime()+(i*minExtensionLength)), Math.max(segment.get5Prime(), spliced_exon.get5Prime()+(i*minExtensionLength)));
						//							System.out.printf("%.2f\t%.2f\t%.2f\t%d\n", extension_read_density/exon_read_density, extension_read_density, exon_read_density, nConsumingReads(bam, strandedness, spliced_exon.chr, s, e, false, segment.isNegativeStrand()));
						written.add(segment.chr, start, end, segment.strand);
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("type", "5p_exon");
						intronic_exon_writer.write("exon", segment.chr, start, end, segment.strand, AnnotatedRegion.GTFAttributeString(attributes));
						//	a					}
					}
				}
			}

			if(closestEnd != null && Math.abs(closestEnd.start - segment.get3Prime()) > minExtensionLength && intervalContains(closestEnd.start, segment.start, segment.end) && !splicedExons.overlappingRegions(segment.chr, segment.get3Prime(), segment.get3Prime(), segment.strand).iterator().hasNext()){
				int s = Math.min(segment.get3Prime(), closestEnd.start-(i*minExtensionLength)); 
				int e = Math.max(segment.get3Prime(), closestEnd.start-(i*minExtensionLength));
				//	a			double extension_read_density = nConsumingReads(bam, strandedness, segment.chr, s, e, false, segment.isNegativeStrand())*1./(e-s+1);
				for(AnnotatedRegion spliced_exon : splicedExons.overlappingRegions(segment.chr, segment.start, segment.end, segment.strand)){
					int start = Math.min(segment.get3Prime(), spliced_exon.get5Prime());
					int end = Math.max(segment.get3Prime(), spliced_exon.get5Prime());
					if(IntervalTools.ContainedIntervals(sj, segment.chr, start, end, segment.strand).size()==0 && !written.contains(segment.chr, start, end, segment.strand)){
						//	a					double exon_read_density = nConsumingReads(bam, strandedness, spliced_exon.chr, spliced_exon.start, spliced_exon.end, false, segment.isNegativeStrand())*1./spliced_exon.size;
						//	a					if(extension_read_density/exon_read_density > minExtensionFraction){
						//							System.out.printf("%s:%d-%d\n", spliced_exon.chr, Math.min(segment.get3Prime(), spliced_exon.get3Prime()-(i*minExtensionLength)), Math.max(segment.get3Prime(), spliced_exon.get3Prime()-(i*minExtensionLength)));
						//							System.out.printf("%.2f\t%.2f\t%.2f\n", extension_read_density/exon_read_density, extension_read_density, exon_read_density);	written.add(segment.chr, s, e, segment.strand);
						written.add(segment.chr, start, end, segment.strand);
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("type", "3p_exon");
						intronic_exon_writer.write("exon", segment.chr, start, end, segment.strand, AnnotatedRegion.GTFAttributeString(attributes));
						//	a					}
					}
				}
			}
		}
	}

	/**
	 * 
	 * @param sfr
	 * @param splice_junction_bed
	 * @param splice_junction_bed2 
	 * @param intronic_exon_gtf 
	 * @param w the distance downstream of the junction to measure coverage
	 * @param d the minimum ratio of spanning:exonic reads that are removed
	 * @param intronic_exon_writer
	 * @throws FileNotFoundException
	 */
	public static void filterIntronicExons(SAMFileReader sfr, Strandedness strandedness, File spliced_exon_gtf, File intronic_exon_gtf, File splice_junction_bed, int w, double d, GTFWriter intronic_exon_writer) throws FileNotFoundException{
		/*
		 * TODO : this method is incomplete, do not use until it is revisited
		 * 
		 */

		BEDIterator sj = new BEDIterator(splice_junction_bed);

		StrandedGenomicIntervalTree<Map<String,Object>> spliced_exons = IntervalTools.buildRegionsTree(new TranscriptIterator(spliced_exon_gtf), true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> intronic_exons = IntervalTools.buildRegionsTree(new TranscriptIterator(intronic_exon_gtf), true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = IntervalTools.buildRegionsTree(sj, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions5p = IntervalTools.buildAttributedTerminiTree(splice_junctions, true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions3p = IntervalTools.buildAttributedTerminiTree(splice_junctions, false, true);

		IntronicExonSupportEvaluater sse = new IntronicExonSupportEvaluater(strandedness, w);

		for(AnnotatedRegion intronic_exon : intronic_exons){
			//			System.out.printf("intronic %s\n", intronic_exon);
			StrandedGenomicIntervalSet intronic_segments = new StrandedGenomicIntervalSet();
			intronic_segments.add(intronic_exon.chr, intronic_exon.start, intronic_exon.end, intronic_exon.strand);
			for(AnnotatedRegion spliced_exon : spliced_exons.overlappingRegions(intronic_exon.chr, intronic_exon.start, intronic_exon.end, intronic_exon.strand)){
				intronic_segments.remove(spliced_exon.chr, spliced_exon.start, spliced_exon.end, spliced_exon.strand);
			}

			for(AnnotatedRegion intronic_segment : intronic_segments){
				if(intronic_segment.size>w){
					//					System.out.println(intronic_exon.toAttributeString());
					//					System.out.printf("seg %s\n", intronic_segment);
					String type = (String) intronic_exon.getAttribute("type");

					if(type.equals("3p_exon")){
						for(AnnotatedRegion splice5p : junctions5p.overlappingRegions(intronic_exon.chr, intronic_exon.start, intronic_exon.end, intronic_exon.strand)){
							sse.reset();
							sse.evaluate(sfr, splice5p.chr, splice5p.start, splice5p.isNegativeStrand(), true);
							splice5p.addAttribute("exon", sse.countExon.n);
							splice5p.addAttribute("intron", sse.countSpanning.n);
							if(sse.countSpanning.n > 30 && sse.countSpanning.n/(0.+sse.countExon.n) > d){
								System.out.printf("%s exon %d intron %d ratio %.2f\n", splice5p, sse.countExon.n, sse.countSpanning.n, sse.countSpanning.n/(0.+sse.countExon.n));
							}
						}
					}
					else if(type.equals("5p_exon")){

					}
				}
			}
		}		

		for(AnnotatedRegion r : junctions5p){
			sse.reset();
			sse.evaluate(sfr, r.chr, r.start, r.isNegativeStrand(), true);
			r.addAttribute("exon", sse.countExon.n);
			r.addAttribute("intron", sse.countSpanning.n);
		}
		for(AnnotatedRegion r : junctions3p){
			sse.reset();
			sse.evaluate(sfr, r.chr, r.start, r.isNegativeStrand(), false);
			r.addAttribute("exon", sse.countExon.n);
			r.addAttribute("intron", sse.countSpanning.n);
		}

		for(AnnotatedRegion r : splice_junctions){
			int nExon5p=-1;
			int nSpan5p=-1;
			int nExon3p=-1;
			int nSpan3p=-1;

			for(AnnotatedRegion j : IntervalTools.OverlappingIntervals(junctions5p, r.chr, r.get5Prime(), r.get5Prime(), r.strand)){
				nExon5p = (Integer) j.getAttribute("exon");
				nSpan5p = (Integer) j.getAttribute("intron");
			}
			for(AnnotatedRegion j : IntervalTools.OverlappingIntervals(junctions3p, r.chr, r.get3Prime(), r.get3Prime(), r.strand)){
				nExon3p = (Integer) j.getAttribute("exon");
				nSpan3p = (Integer) j.getAttribute("intron");
			}

			if(nSpan5p/(0.+nExon5p) < d && nSpan3p/(0.+nExon3p) < d){
				intronic_exon_writer.write(r.chr, r.start, r.end, r.strand);
			}
		}

	}

	public static class IntronicExonSupportEvaluater{

		SAMRecordStrandednessFilter readFilter;
		ConsumingReadFilter consumingExon;
		ConsumingReadFilter consumingIntron;

		Counter<SAMRecord> countExon;
		Counter<SAMRecord> countSpanning;

		int w;

		public IntronicExonSupportEvaluater(Strandedness strandedness, int w) {
			this.w = w;
			readFilter = new SAMRecordStrandednessFilter(strandedness);
			ComposableFilter<SAMRecord> multimap = new MultiMappingfilter(1);
			//			ComposableFilter<SAMRecord> mismatch = new MismatchFilter(0);

			consumingExon = new ConsumingReadFilter(null,0,0);
			consumingIntron = new ConsumingReadFilter(null,0,0);

			countExon = new Counter<SAMRecord>();
			countSpanning = new Counter<SAMRecord>();
			//			Transformer<SAMRecord, String> ss = new Transformer<SAMRecord, String>() {
			//				public String transform(SAMRecord sr) {
			//					List<SAMTagAndValue> l = sr.getAttributes();
			//					for(SAMTagAndValue s : l){
			//						System.out.printf("%s %s\n", s.tag, s.value);
			//					}
			//					return Util.sprintf("%s", sr.getCigarString());
			//				}
			//			};
			//
			//			Echo<String> e = new Echo<String>(System.out);

			readFilter.append(multimap);
			//			readFilter.append(ss);

			//			multimap.append(mismatch);
			//
			//			mismatch.append(consumingExon);

			multimap.append(consumingExon);

			consumingExon.append(countExon);
			consumingExon.append(consumingIntron);

			consumingIntron.append(countSpanning);
			//						consumingIntron.append(ss);

			//						ss.addListener(e);
		}

		public void evaluate(SAMFileReader sfr, String chr, int position, boolean isNegativeStrand, boolean is5p){
			reset();

			int i=is5p^isNegativeStrand ? 1 : 0;
			int j=is5p^isNegativeStrand ? 0 : 1;
			int k=is5p^isNegativeStrand ? -1 : 1;

			readFilter.isNegativeStrand = isNegativeStrand;

			consumingExon.chr = chr;
			consumingExon.start = position + k;
			consumingExon.end = position + k;

			consumingIntron.chr = chr;
			consumingIntron.start = position - (w*k);
			consumingIntron.end = position - (w*k);

			SAMRecordIterator sri = sfr.queryOverlapping(chr,position-i,position+j);
			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				readFilter.add(sr);
			}
			sri.close();		
		}

		private void reset() {
			countExon.reset();
			countSpanning.reset();
		}
	}

	public static void writeUnsplicedRemainder(File segment_bed, File assembly_gtf, GTFWriter remainder_writer) throws FileNotFoundException {
		StrandedGenomicIntervalSet segments = new StrandedGenomicIntervalSet();

		BEDIterator bi = new BEDIterator(segment_bed);
		for(AnnotatedRegion r : bi){
			segments.add(r.chr, r.start, r.end, r.strand);
		}

		TranscriptIterator ti = new TranscriptIterator(assembly_gtf);
		for(AnnotatedRegion r : ti){
			if(r.annotation.equals("exon")){
				segments.remove(r.chr, r.start, r.end, r.strand);
			}
		}

		int segment_count = 0;
		for(AnnotatedRegion segment : segments){
			Map<String,Object> attributes = new HashMap<String, Object>();
			attributes.put("segment_id", Util.sprintf("segment.%08d", segment_count));
			remainder_writer.write("exon", segment.chr, segment.start, segment.end, segment.strand, AnnotatedRegion.GTFAttributeString(attributes));
			segment_count++;
		}
	}

	/*
	 * 	MapList<String, Double[]> comparison = new MapList<String, Double[]>();
		NormalDistribution nd = new NormalDistribution(0, .5);
		NormalDistribution nd2 = new NormalDistribution(3, .5);
		//		for(int n : Util.list(50,100,500,1000,10000,50000,500000,5000000)){
		for(int n : Util.list(50,100,500,1000)){
			// error
			double epsilon = 1e-3;
			// probability that the error is exceeded
			double delta = 1e-2;
			DoubleQuantileFinder qff = QuantileFinderFactory.newDoubleQuantileFinder(true, n, epsilon, delta, 100, RandomEngine.makeDefault());

			for(int i=0; i<n; i++){
				if(Math.random()<.25)
					qff.add(nd.inverseCumulativeProbability(Math.random()));
				else
					qff.add(nd2.inverseCumulativeProbability(Math.random()));
			}

			DoubleArrayList dal = new DoubleArrayList();
			for(int i=0; i<=10; i++){
				dal.add(i*1./10);
			}

			MapList<String, Double[]> data = new MapList<String, Double[]>();

			DoubleArrayList quantiles = qff.quantileElements(dal);
			for (int Ji = 0; Ji < quantiles.size(); Ji++) {
				data.put("list", new Double[]{quantiles.get(Ji),dal.get(Ji)});
			}

			for(double d = -3; d<5; d+=.1){
				double q = qff.phi(d);
				data.put("exact", new Double[]{d,(.25*nd.cumulativeProbability(d))+(.75*nd2.cumulativeProbability(d))});

				double e = (.25*nd.cumulativeProbability(d))+(.75*nd2.cumulativeProbability(d));

				data.put("estimate", new Double[]{d,q});

				//			System.out.printf("q %.2f e %.2f\n", q,e);
				comparison.put(""+n, new Double[]{q,e-q});
			}

			//		System.out.println(qff.memory());
			//		System.out.println(qff.totalMemory());
			ChartUtils.showChart(ChartUtils.createScatterChart(data.getMap(), "", "x", "cdf", false));
		}

	 */

	/*
	 * boolean known_upper_bound = true;
		int nSamples = 10000;
		NormalDistribution nd = new NormalDistribution(0, 1);
		double nsd = .5;
		double d1 = nd.cumulativeProbability(nsd);
		double d2 = nd.cumulativeProbability(-nsd);
		System.out.println(d1);
		NormalDistribution nd1 = new NormalDistribution(10, 8);
		NormalDistribution nd2 = new NormalDistribution(1000, 300);
		DoubleQuantileFinder dqf= QuantileFinderFactory.newDoubleQuantileFinder(known_upper_bound, nSamples, 1e-3, 1e-3, 1000, RandomEngine.makeDefault());
		MapList<String, Double> dat = new MapList<String,Double>();
		MomentTracker mt = new MomentTracker();
		for(int i=0; i<nSamples; i++){
			Double db = null;
			if(Math.random()<.99){
				db = nd1.inverseCumulativeProbability(Math.random());
			}
			else{
				db = nd2.inverseCumulativeProbability(Math.random());
			}
			mt.add(db);
			dqf.add(db);
			dat.put("", db);
		}

		DoubleArrayList dal = new DoubleArrayList(2);
		dal.add(d2);
		dal.add(.5);
		dal.add(d1);

		DoubleArrayList q = dqf.quantileElements(dal);
		System.out.println(q.get(0));
		System.out.println(q.get(1));

		double mean = (q.get(1)+q.get(2)+q.get(0))/3;
		double sd1 = (q.get(2)-mean)/nsd;
		double sd2 = (mean-q.get(0))/nsd;

	 */

	public static int estimateMateInsertSize(SAMFileReader sfr, Strandedness strandedness, File segment_bed, double quantile) throws FileNotFoundException{
		BEDIterator ti = new BEDIterator(segment_bed);

		NormalDistribution nd = new NormalDistribution(0, 1);
		// number of standard deviations to use quantile to estimate robust standard deviation 
		// effective if ~0.841 percent of reads map at proper distance [this is P(N(0,1) < 1.0)]
		double nsd = 1;

		double dneg = nd.cumulativeProbability(-nsd);
		double dzer = nd.cumulativeProbability(0);
		double dpos = nd.cumulativeProbability(nsd);

		//		List<Integer> insert_sizes = new LinkedList<Integer>();
		StrandedGenomicIntervalSet intervals = IntervalTools.buildStrandedIntervalSet(ti);

		boolean known_upper_bound = true;
		DoubleQuantileFinder dqf= QuantileFinderFactory.newDoubleQuantileFinder(known_upper_bound, BAMTools.totalAlignedReads(sfr), 1e-2, 1e-2, 1000, RandomEngine.makeDefault());

		SAMRecordIterator sri = sfr.iterator();

		while(sri.hasNext()){
			SAMRecord sr = sri.next();

			String chr = sr.getReferenceName();

			boolean mateUnmapped = sr.getMateUnmappedFlag();
			boolean isNegativeStrand = sr.getReadNegativeStrandFlag();
			int nHits = sr.getIntegerAttribute("NH");
			if(!isNegativeStrand && !mateUnmapped && nHits==1){
				String mateChr = sr.getMateReferenceName();
				int readEnd = sr.getAlignmentEnd();
				int mateStart = sr.getMateAlignmentStart();

				// if the mate is contained in the same continuous segment
				if(chr.equals(mateChr) && IntervalTools.isContained(intervals, chr, Math.min(readEnd,mateStart), Math.max(readEnd,mateStart), BAMTools.strand(sr, strandedness))){
					int insert_size = mateStart-readEnd-1;

					//						insert_sizes.add(insert_size);
					dqf.add(insert_size);
				}
			}
			//				if(insert_sizes.size()>1000000)
			//					break;
		}

		sri.close();


		DoubleArrayList dal = new DoubleArrayList(3);
		dal.add(dneg);
		dal.add(dzer);
		dal.add(dpos);

		DoubleArrayList q = dqf.quantileElements(dal);


		double mean = q.get(1);
		// estimate sd from the average of two quantiles
		double sd1 = (q.get(2)-mean)/nsd;
		double sd2 = (mean-q.get(0))/nsd;
		NormalDistribution estimated_insert_distribution = new NormalDistribution(mean, (sd1+sd2)/2);
		return (int) estimated_insert_distribution.inverseCumulativeProbability(quantile);
	}



	public static void identifySpannableRegions(SAMFileReader sfr, Strandedness strandedness, BEDWriter matepair_bw, int max_insert_size) throws FileNotFoundException{
		StrandedGenomicIntervalSet mate_spanned_regions = new StrandedGenomicIntervalSet();

		SAMRecordIterator sri = sfr.iterator();
		//		SAMRecordIterator sri = sfr.queryOverlapping("19", 3892522, 3898022);
		//sfr.iterator();
		while(sri.hasNext()){

			SAMRecord sr = sri.next();
			//			System.out.println(sr);
			String chr = sr.getReferenceName();
			int end = sr.getAlignmentEnd();

			boolean unmapped = sr.getReadUnmappedFlag();
			boolean mateUnmapped = sr.getMateUnmappedFlag();
			boolean isNegativeStrand = sr.getReadNegativeStrandFlag();

			char strand = BAMTools.strand(sr, strandedness);
			if(!unmapped && !isNegativeStrand && !mateUnmapped){
				//				System.out.println("there");
				String mate_chr = sr.getMateReferenceName();
				int mate_start = sr.getMateAlignmentStart();
				int insert_size = mate_start-sr.getAlignmentEnd()-1;
				//				System.out.println(insert_size);
				if(chr.equals(mate_chr) && insert_size>0 && insert_size <= max_insert_size){
					mate_spanned_regions.add(sr.getReferenceName(), end+1, mate_start-1, strand);
				}
			}
		}

		sri.close();

		for(AnnotatedRegion r : mate_spanned_regions){
			matepair_bw.write(r.chr, r.start, r.end, r.strand);
		}

	}

	public static void scaffoldSpannableRegions(File segment_bed,	File matepair_bed, BEDWriter scaffolded_bw) throws FileNotFoundException {
		BEDIterator bi = new BEDIterator(segment_bed);
		StrandedGenomicIntervalSet t = IntervalTools.buildStrandedIntervalSet(bi);
		IntervalTools.addRegions(t, new BEDIterator(matepair_bed));

		for(AnnotatedRegion r : t){
			scaffolded_bw.write(r.chr, r.start, r.end, r.strand);
		}
	}

	public static void main(String[] args) throws FileNotFoundException {
		if(true){
			SAMFileReader sfr = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE51572/test/g.bam"));
			File segment_bed = new File("/mnt/LaiLab/sol/GSE51572/test/isoscm/tmp/gu.seg.bed");
//			countJunctionSupportingReads(sfr, new File("/mnt/LaiLab/sol/GSE51572/test/isoscm/tmp/g.sj.bed"), new GTFWriter("/mnt/LaiLab/sol/GSE51572/test/g1.counts.gtf"));

			int i = estimateMateInsertSize(sfr, Strandedness.unstranded, segment_bed, .99);
			System.out.println(i);
		}
		if(false)
		{
			SAMFileReader sfr = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594413.bam"));
			BEDWriter matepair_bw = new BEDWriter("/mnt/LaiLab/sol/GSE41637/mapped/test_pairs.bed");
			identifySpannableRegions(sfr, Strandedness.reverse_forward, matepair_bw, 160);
			matepair_bw.close();
		}
		if(false)
		{
			trimInferredExons(new File("/home/sol/lailab/sol/sangercenter/de_dup/isoscm/gtf/tmp/hip.exon.gtf"), new GTFWriter("/home/sol/lailab/sol/sangercenter/de_dup/isoscm/gtf/tmp/hip.exon.trimmed.gtf"));
		}
		if(false)
		{
			BEDWriter bw = new BEDWriter("/mnt/LaiLab/sol/sangercenter/de_dup/isoscm/gtf/tmp/hip.inferred.bed");
			inferSegmentStrand(new File("/mnt/LaiLab/sol/sangercenter/de_dup/isoscm/gtf/tmp/hip.sj.bed"), new File("/home/sol/workspace/IsoSCM/unstranded.bed"),bw);
			bw.close();
		}
		if(false)
		{			
			File segment_bed = new File("/mnt/LaiLab/sol/GSE41637/nb/gtf/tmp/SRR594393.seg.bed");
			File repeat_bed = new File("/home/sol/data/repeats/mm10/mm10.RepeatMasker.bed");
			BEDWriter merged_segments = new BEDWriter("/mnt/LaiLab/sol/GSE41637/nb/test.merge.repeat.bed");
			mergeSegments(segment_bed, repeat_bed, merged_segments, 100);
			merged_segments.close();
		}
		if(false)
		{
			File segment_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.seg.bed");
			BEDWriter merged_segments = new BEDWriter("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.seg.merged.bed");
			//			mergeSegments(segment_bed, merged_segments, 100);
			merged_segments.close();
		}
		if(false)
		{
			File splice_junction_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.sj.filtered.bed");
			File spliced_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.exon.gtf");
			File intronic_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/test.intronic.gtf");
			SAMFileReader bam = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			Strandedness strandedness = Strandedness.reverse_forward;
			int w = 5;
			double d = .1;

			GTFWriter intronic_exon_writer = new GTFWriter("/home/sol/workspace/IsoSCM/tests/gtf/tmp/test.intronic.filtered.gtf");
			filterIntronicExons(bam, strandedness, spliced_exon_gtf, intronic_exon_gtf, splice_junction_bed, w, d, intronic_exon_writer);
			intronic_exon_writer.close();
		}
		if(false)
		{
			File segment_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.seg.bed");
			File splice_junction_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.sj.bed");
			File spliced_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.exon.gtf");
			SAMFileReader bam = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			Strandedness strandedness = Strandedness.reverse_forward;
			int minExtensionLength = 5;
			double minExtensionFraction = .1;
			GTFWriter intronic_exon_writer = new GTFWriter("/home/sol/workspace/IsoSCM/tests/gtf/tmp/test.intronic.gtf");
			identifyIntronicExons(segment_bed, splice_junction_bed, spliced_exon_gtf, bam, strandedness, intronic_exon_writer, minExtensionLength, minExtensionFraction);
			intronic_exon_writer.close();
		}
		System.exit(0);

		//				GTFWriter gw = new GTFWriter("test5.consolidated.gtf");
		//		GTFWriter gw = new GTFWriter("test6.consolidated.gtf");
		//				intervals("test5.gtf",gw);
		SAMFileReader sfr = new SAMFileReader(new File("/home/sol/data/sangercenter/hippocampus.bam"));

		ConsolidatedIntervalResult expressedSegmentsResult = expressedSegmentIntervals("test6.gtf");
		GenomicIntervalTree<Map<String,Object>> expressedSegments = expressedSegmentsResult.segments;
		GenomicIntervalSet highlyExpressedSegments = new GenomicIntervalSet();
		double expressed_threshold = 10;
		for(AnnotatedRegion r : expressedSegments){
			if((Double)r.getAttribute("mle")>expressed_threshold){
				highlyExpressedSegments.add(r.chr, r.start, r.end);
			}
		}

		System.out.println(expressedSegments.size());
		System.out.println(highlyExpressedSegments.size());


		StrandedGenomicIntervalTree<Map<String,Object>> sj = new StrandedGenomicIntervalTree<Map<String,Object>>();
		BEDIterator splice_junctions = new BEDIterator("hippocampus.sj.bed");
		for(AnnotatedRegion splice_junction : splice_junctions){
			//			System.out.println(splice_junction);
			if(!sj.contains(splice_junction.chr, splice_junction.start, splice_junction.end, splice_junction.strand))
				sj.add(splice_junction);
		}

		int w = 50;
		int hw = (w+1)/2;

		StrandedGenomicIntervalSet clusteredSegments = new StrandedGenomicIntervalSet();
		for(AnnotatedRegion sje : sj){
			boolean overlaps = false;
			for(AnnotatedRegion r : highlyExpressedSegments.overlappingRegions(sje.chr, sje.start-w-hw, sje.start+w-hw)){
				overlaps |= true;
			}
			for(AnnotatedRegion r : highlyExpressedSegments.overlappingRegions(sje.chr, sje.end-w-hw, sje.end+w-hw)){
				overlaps |= true;
			}
			if(overlaps){
				clusteredSegments.add(sje.chr, sje.start-w, sje.end+w, sje.strand);
			}
		}

		/**
		 * TODO rather than have unstranded clustered segments, make them stranded above. Then, only add highly expressed segments that
		 * overlap with these.
		 */

		for(AnnotatedRegion r : highlyExpressedSegments){
			boolean overlaps_splice = false;
			for(char strand : Util.list('+','-')){
				if(clusteredSegments.overlappingRegions(r.chr, r.start, r.end, strand).iterator().hasNext()){
					clusteredSegments.add(r.chr, r.start, r.end, strand);
					overlaps_splice |= true;
				}
			}
			if(!overlaps_splice)
				clusteredSegments.add(r.chr, r.start, r.end, '.');

		}
		//		chr1:23367501-23929101
		BEDWriter bw = new BEDWriter("test6.clusters.bed");
		//		BEDWriter models = new BEDWriter("test6.models.bed");
		BEDWriter plus_models = new BEDWriter("test6.plus_models.bed");
		BEDWriter minus_models = new BEDWriter("test6.minus_models.bed");
		BEDWriter bw_i = new BEDWriter("increase.bed");
		BEDWriter bw_d = new BEDWriter("decrease.bed");
		GTFWriter transcribed_fragments = new GTFWriter("test6.transfrags.gtf");
		GTFWriter assembly = new GTFWriter("chr1.assembly.gtf");

		int cluster_id=0;
		for(AnnotatedRegion r : clusteredSegments){
			cluster_id++;
			bw.write(r.chr, r.start, r.end, r.strand);

			// save the first and last increase and decrease, these should not be removed even if they are near a splice junction
			ExtremeObjectTracker<AnnotatedRegion, Integer> terminalIncreaseTracker = new ExtremeObjectTracker<AnnotatedRegion, Integer>(new Util.ComparableComparator<Integer>());
			ExtremeObjectTracker<AnnotatedRegion, Integer> terminalDecreaseTracker = new ExtremeObjectTracker<AnnotatedRegion, Integer>(new Util.ComparableComparator<Integer>());

			GenomicIntervalTree<Map<String,Object>> increasePoints = new GenomicIntervalTree<Map<String,Object>>();
			GenomicIntervalTree<Map<String,Object>> decreasePoints = new GenomicIntervalTree<Map<String,Object>>();

			for(AnnotatedRegion changepoint : expressedSegmentsResult.increasePoints.overlappingRegions(r.chr, r.start, r.end)){
				terminalIncreaseTracker.put(changepoint, changepoint.start);
				increasePoints.add(changepoint.chr, changepoint.start, changepoint.end, changepoint.attributes);
			}
			for(AnnotatedRegion changepoint : expressedSegmentsResult.decreasePoints.overlappingRegions(r.chr, r.start, r.end)){
				terminalDecreaseTracker.put(changepoint, changepoint.start);
				decreasePoints.add(changepoint.chr, changepoint.start, changepoint.end, changepoint.attributes);
			}

			AnnotatedRegion terminalIncrease = terminalIncreaseTracker.getMinObject();
			AnnotatedRegion terminalDecrease = terminalDecreaseTracker.getMaxObject();

			List<AnnotatedRegion> to_remove = new LinkedList<AnnotatedRegion>();
			for(char strand : Util.list('+','-')){
				for(AnnotatedRegion sje : sj.overlappingRegions(r.chr, r.start, r.end, strand)){
					for(AnnotatedRegion changepoint : expressedSegmentsResult.changePoints.overlappingRegions(sje.chr, sje.start-w-hw, sje.start+w-hw)){
						to_remove.add(changepoint);
					}
					for(AnnotatedRegion changepoint : expressedSegmentsResult.changePoints.overlappingRegions(sje.chr, sje.end-w-hw, sje.end+w-hw)){
						to_remove.add(changepoint);
					}
				}
			}

			for(AnnotatedRegion changepoint : to_remove){
				if( (terminalIncrease!=null && changepoint.start != terminalIncrease.start) && (terminalDecrease!=null && changepoint.start != terminalDecrease.start)){
					increasePoints.remove(changepoint.chr, changepoint.start, changepoint.end);
					decreasePoints.remove(changepoint.chr, changepoint.start, changepoint.end);
				}
			}

			List<AnnotatedRegion> increases = new LinkedList<AnnotatedRegion>();
			List<AnnotatedRegion> decreases = new LinkedList<AnnotatedRegion>();

			for(AnnotatedRegion changepoint : increasePoints.overlappingRegions(r.chr, r.start, r.end)){
				if( (Double) changepoint.getAttribute("after_mle") > expressed_threshold ){
					increases.add(changepoint);
					bw_i.write(changepoint.chr, changepoint.start, changepoint.end, changepoint.strand);
				}
			}
			for(AnnotatedRegion changepoint : decreasePoints.overlappingRegions(r.chr, r.start, r.end)){
				if( (Double) changepoint.getAttribute("before_mle") > expressed_threshold ){
					decreases.add(changepoint);
					bw_d.write(changepoint.chr, changepoint.start, changepoint.end, changepoint.strand);
				}
			}

			//			if(r.strand != '.')
			//				System.out.println(r);
			//			for(AnnotatedRegion increase : increases){
			//				for(AnnotatedRegion decrease : decreases){
			//					if(increase.start < decrease.end){
			//
			//						List<AnnotatedRegion> splices = prepareSpliceJunctions(sfr, sj, r.chr, increase.start, decrease.end, r.strand);
			//
			//						enumerateExons(r.chr, increase.start, decrease.end, r.strand, splices, models);
			//
			//						//					bw.write(r.chr, increase.start, decrease.end, '.');
			//					}
			//				}
			//			}

			//TODO make GTF output with junctions, exons 5' and 3' labeled

			if(r.strand != '.'){
				BEDWriter models = r.isNegativeStrand()? minus_models : plus_models;

				System.out.println(r);

				StrandedGenomicIntervalTree<Map<String,Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();

				StrandedGenomicPositionTree splice5ps = new StrandedGenomicPositionTree();
				StrandedGenomicPositionTree splice3ps = new StrandedGenomicPositionTree();

				// find all splice junctions
				//				List<AnnotatedRegion> splices = prepareSpliceJunctions(sfr, sj, r.chr, r.start, r.end, r.strand);
				for(AnnotatedRegion splice : sj.overlappingRegions(r.chr, r.start, r.end, r.strand)){
					//				for(AnnotatedRegion splice : splices){
					StrandedGenomicIntervalSet exonIntervals = new StrandedGenomicIntervalSet();
					exonIntervals.add(splice.chr, splice.start-1, splice.start-1, splice.strand);
					exonIntervals.add(splice.chr, splice.end+1, splice.end+1, splice.strand);
					writeTranscript(exonIntervals, models);
					//					models.write(splice.chr, splice.start, splice.end, splice.strand);
					splice5ps.add(splice.chr, splice.get5Prime(), splice.strand);
					splice3ps.add(splice.chr, splice.get3Prime(), splice.strand);

					Map<String,Object> attributes = new HashMap<String, Object>();
					attributes.put("gene_id", splice.toString());
					attributes.put("transcript_id", splice.toString());
					attributes.put("type", "splice_junction");

					assembly.write("exon", splice.chr, splice.start-1, splice.start-1, splice.strand, AnnotatedRegion.GTFAttributeString(attributes));
					assembly.write("exon", splice.chr, splice.end+1, splice.end+1, splice.strand, AnnotatedRegion.GTFAttributeString(attributes));
				}

				// find all internal exons
				for(AnnotatedRegion splice : sj.overlappingRegions(r.chr, r.start, r.end, r.strand)){
					AnnotatedRegion exon5p = splice3ps.getClosestUpstream(splice.chr, splice.get5Prime(), splice.strand);
					if(exon5p!=null){
						int s5p = splice.get5Prime();
						int e5p = exon5p.start;
						if(!exons.contains(r.chr, Math.min(s5p,e5p)+1, Math.max(s5p,e5p)-1, r.strand)){
							{
								Map attributes = new HashMap();
								attributes.put("type", "internal");
								transcribed_fragments.write("exon", r.chr, Math.min(s5p,e5p)+1, Math.max(s5p,e5p)-1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
								assembly.write("exon", r.chr, Math.min(s5p,e5p)+1, Math.max(s5p,e5p)-1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}

							models.write(r.chr, Math.min(s5p,e5p)+1, Math.max(s5p,e5p)-1, r.strand);
							exons.add(r.chr, Math.min(s5p,e5p)+1, Math.max(s5p,e5p)-1, r.strand);
						}
					}
				}

				// find all terminal exons
				List<AnnotatedRegion> ends3p = r.isNegativeStrand()?increases:decreases;
				List<AnnotatedRegion> ends5p = !r.isNegativeStrand()?increases:decreases;

				// the 3' exon
				// since there will be some error due to granularity of SCM, don't forget to add a half-width
				// splice junctions also have to be adjusted by one
				int i3p = r.isNegativeStrand()?-1:1;
				int i5p = !r.isNegativeStrand()?-1:1;

				for(AnnotatedRegion end3p : ends3p){

					AnnotatedRegion transcriptBoundary = splice3ps.getClosestUpstream(r.chr, end3p.start, r.strand);
					if(transcriptBoundary!=null){
						int s5p = end3p.start;
						int e5p = transcriptBoundary.start;
						if(!exons.contains(r.chr, Math.min(s5p,e5p), Math.max(s5p,e5p), r.strand)){

							{	
								Map attributes = new HashMap();
								attributes.put("type", "exon3p");
								transcribed_fragments.write("exon", r.chr, Math.min(s5p+hw,e5p)+1, Math.max(s5p+hw,e5p)-1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
							{
								Map attributes = new HashMap();
								attributes.put("type", "3'end");
								assembly.write("exon", r.chr, s5p+hw, s5p+hw, r.strand, AnnotatedRegion.GTFAttributeString(attributes));	
							}

							models.write(r.chr, Math.min(s5p+hw,e5p+i3p), Math.max(s5p+hw,e5p+i3p), r.strand);
							exons.add(r.chr, Math.min(s5p,e5p), Math.max(s5p,e5p), r.strand);
						}
					}
				}
				// the 5' exon
				for(AnnotatedRegion end5p : ends5p){

					AnnotatedRegion transcriptBoundary = splice5ps.getClosestDownstream(r.chr, end5p.start, r.strand);
					if(transcriptBoundary!=null){
						int s5p = end5p.start;
						int e5p = transcriptBoundary.start;
						if(!exons.contains(r.chr, Math.min(s5p,e5p), Math.max(s5p,e5p), r.strand)){
							{
								Map attributes = new HashMap();
								attributes.put("type", "exon5p");
								transcribed_fragments.write("exon", r.chr, Math.min(s5p+hw,e5p)+1, Math.max(s5p+hw,e5p)-1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
							{
								Map attributes = new HashMap();
								attributes.put("type", "5'end");
								assembly.write("exon", r.chr, s5p+hw, s5p+hw, r.strand, AnnotatedRegion.GTFAttributeString(attributes));	
							}
							models.write(r.chr, Math.min(s5p+hw,e5p+i5p), Math.max(s5p+hw,e5p+i5p), r.strand);
							exons.add(r.chr, Math.min(s5p,e5p), Math.max(s5p,e5p), r.strand);
						}
					}
				}

				//				for(AnnotatedRegion increase : increases){
				//					for(AnnotatedRegion decrease : decreases){
				//						if(increase.start < decrease.end){
				//
				//							List<AnnotatedRegion> splices = prepareSpliceJunctions(sfr, sj, r.chr, increase.start, decrease.end, r.strand);
				//
				//							enumerateExons(r.chr, increase.start, decrease.end, r.strand, splices, models);
				//
				//							//					bw.write(r.chr, increase.start, decrease.end, '.');
				//						}
				//					}
				//				}
			}
		}
		bw_i.close();
		bw_d.close();
		bw.close();
		transcribed_fragments.close();
		//					gw.close();
	}


}
