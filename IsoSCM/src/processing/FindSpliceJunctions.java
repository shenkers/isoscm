package processing;

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Array;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Map.Entry;
import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeMap;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.math3.special.Beta;

import tools.AnnotatedRegion;
import tools.BAMTools;
import tools.BEDTools.BEDWriter;
import tools.GTFTools.GTFWriter;
import tools.IntervalTools;
import tools.ParseBed.BEDIterator;
import tools.ParseGTF.TranscriptIterator;
import tools.StrandedGenomicIntervalSet;
import tools.StrandedGenomicIntervalTree;
import tools.Strandedness;
import util.Util;
import util.Util.ExtremeObjectTracker;
import util.Util.Instantiator;
import util.Util.MapCounter;
import util.Util.MapFactory;
import filter.ComposableFilter;
import filter.ConsumingReadFilter;
import filter.Counter;
import filter.MultiMappingfilter;
import filter.SAMRecordStrandednessFilter;

public class FindSpliceJunctions {

	public static List<AnnotatedRegion> spliceJunction(SAMFileReader sfr, String chr, int start, int end, boolean contained){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		List<AnnotatedRegion> splice_junctions = new LinkedList<AnnotatedRegion>();

		while(sri.hasNext()){
			SAMRecord sr = sri.next();
			if(sr.getCigarString().indexOf('N') == -1)
				continue;

			Cigar cigar = sr.getCigar();
			int alignmentPosition = sr.getAlignmentStart();

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
					if(isIntron){
						//							System.out.println(Util.list(alignmentPosition,alignmentPosition+cigarElement.getLength()-1));
						if((alignmentPosition >= start && alignmentPosition <= end) || (alignmentPosition+cigarElement.getLength()-1 >= start && alignmentPosition+cigarElement.getLength()-1 <= end))
							splice_junctions.add(new AnnotatedRegion("splice_junction", sr.getReferenceName(), alignmentPosition,alignmentPosition+cigarElement.getLength()-1, sr.getCharacterAttribute("XS")));
						//							System.out.println(alignmentPosition);
					}
					alignmentPosition += cigarElement.getLength();
					//						if(isIntron
					////								&& alignmentPosition >= start && alignmentPosition <= end
					//								){
					//							n++;
					//							System.out.println(alignmentPosition);
					//							System.out.println();
					//							break spliced;
					//						}
				}			
			}

		}
		sri.close();

		return splice_junctions;
	}

	public static List<AnnotatedRegion> splice_junctions(SAMRecord sr){
		List<AnnotatedRegion> splice_junctions = new LinkedList<AnnotatedRegion>();

		if(sr.getCigarString().indexOf('N') != -1){

			Cigar cigar = sr.getCigar();
			int alignmentPosition = sr.getAlignmentStart();

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
					if(isIntron){
						AnnotatedRegion sj = new AnnotatedRegion("sj", sr.getReferenceName(), alignmentPosition,alignmentPosition+cigarElement.getLength()-1, sr.getCharacterAttribute("XS"));
						splice_junctions.add(sj);
					}
					alignmentPosition += cigarElement.getLength();
				}			
			}

		}

		return splice_junctions;
	}

	public static StrandedGenomicIntervalTree<Map<String,Object>> tabulateSpliceJunctions(SAMFileReader sfr, Strandedness strandedness, BEDWriter bw){
		Map<String,Integer> referenceLengths = BAMTools.referenceSequenceLengths(sfr.getFileHeader());

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(String chr : referenceLengths.keySet()){
			int start = 1;
			int end = referenceLengths.get(chr);
			boolean contained = true;

			SAMRecordIterator sri = sfr.query(chr, start, end, contained);
			//		List<AnnotatedRegion> splice_junctions = new LinkedList<AnnotatedRegion>();
			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				if(sr.getCigarString().indexOf('N') == -1)
					continue;
				
				if(strandedness!=Strandedness.unstranded && sr.getAttribute("XS")!=null && BAMTools.strand(sr, strandedness)!=sr.getCharacterAttribute("XS"))
					continue;

				Cigar cigar = sr.getCigar();
				int alignmentPosition = sr.getAlignmentStart();

				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
						if(isIntron){
							//							System.out.println(Util.list(alignmentPosition,alignmentPosition+cigarElement.getLength()-1));
							//							if((alignmentPosition >= start && alignmentPosition <= end) || (alignmentPosition+cigarElement.getLength()-1 >= start && alignmentPosition+cigarElement.getLength()-1 <= end)){
							AnnotatedRegion sj = new AnnotatedRegion("sj", chr, alignmentPosition,alignmentPosition+cigarElement.getLength()-1, sr.getCharacterAttribute("XS"));
							if(!splice_junctions.contains(sj.chr, sj.start, sj.end, sj.strand)){
								splice_junctions.add(sj);
								bw.write(sj.chr, sj.start, sj.end, sj.strand);
							}
							//							}
							//							System.out.println(alignmentPosition);
						}
						alignmentPosition += cigarElement.getLength();
						//						if(isIntron
						////								&& alignmentPosition >= start && alignmentPosition <= end
						//								){
						//							n++;
						//							System.out.println(alignmentPosition);
						//							System.out.println();
						//							break spliced;
						//						}
					}			
				}

			}
			sri.close();
		}

		return splice_junctions;
	}

	/**
	 * 
	 * @param sfr
	 * @param splice_junction_bed
	 * @param w the distance downstream of the junction to measure coverage
	 * @param d the minimum ratio of spanning:exonic reads that are removed
	 * @param bw
	 * @throws FileNotFoundException
	 */
	public static void filterSpliceJunctions(SAMFileReader sfr, Strandedness strandedness, File splice_junction_bed, int w, double d, BEDWriter bw) throws FileNotFoundException{
		BEDIterator sj = new BEDIterator(splice_junction_bed);

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = IntervalTools.buildRegionsTree(sj, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions5p = IntervalTools.buildAttributedTerminiTree(splice_junctions, true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions3p = IntervalTools.buildAttributedTerminiTree(splice_junctions, false, true);

		SpliceSupportEvaluater sse = new SpliceSupportEvaluater(strandedness, w);

		for(AnnotatedRegion r : junctions5p){
			sse.reset();
			sse.evaluate(sfr, r.chr, r.start, r.isNegativeStrand(), true);
			r.addAttribute("exon", sse.countExon.n);
			r.addAttribute("intron", sse.countSpanning.n);

			//			System.out.println(r);
			//			System.out.println(r.toAttributeString());

			//			int nExon = sse.countExon.n;
			//			int nSpanning = sse.countSpanning.n;
			//			
			//			if(nExon+nSpanning!=0 && nSpanning/(0.+nExon) > .4 && nSpanning/(0.+nExon) <= .5){
			////				System.out.printf("retained\n", r);
			//				System.out.printf("%s %.2f %d %d\n", r, nSpanning/(0.+nExon), nExon, nSpanning);
			//			}
			//			if(nExon+nSpanning!=0){
			////				System.out.printf("retained\n", r);
			////				System.out.printf("%s %.2f %d %d\n", r, nSpanning/(0.+nExon+nSpanning), nExon, nSpanning);
			//			}
		}
		for(AnnotatedRegion r : junctions3p){
			sse.reset();
			sse.evaluate(sfr, r.chr, r.start, r.isNegativeStrand(), false);
			r.addAttribute("exon", sse.countExon.n);
			r.addAttribute("intron", sse.countSpanning.n);
		}

		//		for(AnnotatedRegion r : junctions5p){
		//			System.out.printf("5p %s %s\n", r, r.toAttributeString());
		//		}
		//		for(AnnotatedRegion r : junctions3p){
		//			System.out.printf("3p %s %s\n", r, r.toAttributeString());
		//		}

		//		StrandedGenomicIntervalTree<Map<String,Object>> satisfy_both = new StrandedGenomicIntervalTree<Map<String,Object>>();
		//		StrandedGenomicIntervalTree<Map<String,Object>> satisfy_5p = new StrandedGenomicIntervalTree<Map<String,Object>>();
		//		StrandedGenomicIntervalTree<Map<String,Object>> satisfy_3p = new StrandedGenomicIntervalTree<Map<String,Object>>();

		// minimum number of reads to perform filtering
		int min_reads = 16;

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

			//			System.out.printf("%s %.2f %.2f %d %d %d %d\n", r, nSpan5p/(0.+nExon5p), nSpan3p/(0.+nExon3p), nSpan5p,nExon5p, nSpan3p,nExon3p);

			//			if((nSpan5p+nExon5p <= min_reads || nSpan5p/(0.+nExon5p) <= d) && (nSpan5p+nExon5p <= min_reads || nSpan3p/(0.+nExon3p) <= d)){
			if((nSpan5p+nExon5p <= min_reads || nSpan5p/(0.+nExon5p) <= d) && (nSpan3p+nExon3p <= min_reads || nSpan3p/(0.+nExon3p) <= d)){
				bw.write(r.chr, r.start, r.end, r.strand);
				//				satisfy_both.add(r.chr,r.start,r.end,r.strand);
			}
			//			else{
			//				if(nSpan5p/(0.+nExon5p) < d){
			//					satisfy_5p.add(r.chr,r.start,r.end,r.strand);
			//				}
			//				if(nSpan3p/(0.+nExon3p) < d){
			//					satisfy_3p.add(r.chr,r.start,r.end,r.strand);
			//				}
			//			}
		}

	}

	public static class SpliceSupportEvaluater{

		SAMRecordStrandednessFilter readFilter;
		ConsumingReadFilter consumingExon;
		ConsumingReadFilter consumingIntron;

		Counter<SAMRecord> countExon;
		Counter<SAMRecord> countSpanning;

		int w;

		public SpliceSupportEvaluater(Strandedness strandedness, int w) {
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

			//			readFilter.append(multimap);
			//			multimap.append(consumingExon);

			readFilter.append(consumingExon);

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

	public static class tSpliceSupportEvaluater{

		SAMRecordStrandednessFilter readFilter;
		ConsumingReadFilter consumingExonW;
		ConsumingReadFilter consumingIntronW;
		ConsumingReadFilter consumingExon;
		ConsumingReadFilter consumingIntron;

		Counter<SAMRecord> countExonW;
		Counter<SAMRecord> countIntronW;
		Counter<SAMRecord> countExon;
		Counter<SAMRecord> countIntron;

		int w;

		public tSpliceSupportEvaluater(Strandedness strandedness, int w) {
			this.w = w;
			readFilter = new SAMRecordStrandednessFilter(strandedness);
			ComposableFilter<SAMRecord> multimap = new MultiMappingfilter(1);

			consumingExonW = new ConsumingReadFilter(null,0,0);
			consumingIntronW = new ConsumingReadFilter(null,0,0);
			consumingExon = new ConsumingReadFilter(null,0,0);
			consumingIntron = new ConsumingReadFilter(null,0,0);

			countExonW = new Counter<SAMRecord>();
			countIntronW = new Counter<SAMRecord>();
			countExon = new Counter<SAMRecord>();
			countIntron = new Counter<SAMRecord>();

			readFilter.append(multimap);

			readFilter.append(consumingExonW);
			readFilter.append(consumingIntronW);
			readFilter.append(consumingExon);
			readFilter.append(consumingIntron);

			consumingExonW.append(countExonW);
			consumingIntronW.append(countIntronW);
			consumingExon.append(countExon);
			consumingIntron.append(countIntron);
		}

		public void evaluate(SAMFileReader sfr, String chr, int position, boolean isNegativeStrand, boolean is5p){
			reset();

			int[] exon = is5p ? IntervalTools.offsetInterval(position, w, -1, isNegativeStrand) : IntervalTools.offsetInterval(position, -1, w, isNegativeStrand);
			//			int[] intron = !is5p ? IntervalTools.offsetInterval(position, 1, 0, isNegativeStrand) : IntervalTools.offsetInterval(position, 0, w, isNegativeStrand);

			readFilter.isNegativeStrand = isNegativeStrand;

			consumingExonW.chr = chr;
			consumingExonW.start = exon[0];
			consumingExonW.end = exon[1];

			consumingIntronW.chr = chr;
			consumingIntronW.start = position;
			consumingIntronW.end = position;

			consumingExon.chr = chr;
			consumingExon.start = IntervalTools.offsetPosition(position, 1, isNegativeStrand, is5p);
			consumingExon.end = IntervalTools.offsetPosition(position, 1, isNegativeStrand, is5p);

			consumingIntron.chr = chr;
			consumingIntron.start = position;
			consumingIntron.end = position;

			SAMRecordIterator sri = sfr.queryOverlapping(chr,position-w,position+w);
			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				readFilter.add(sr);
			}
			sri.close();		
		}

		private void reset() {
			countExonW.reset();
			countIntronW.reset();
			countExon.reset();
			countIntron.reset();
		}
	}

	public static void tFilterSpliceJunctions(SAMFileReader sfr, Strandedness strandedness, File splice_junction_bed, int w, double alpha, GTFWriter gw) throws FileNotFoundException{
		BEDIterator sj = new BEDIterator(splice_junction_bed);

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = IntervalTools.buildRegionsTree(sj, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> filtered_junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();//.buildRegionsTree(sj, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions5p = IntervalTools.buildAttributedTerminiTree(splice_junctions, true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions3p = IntervalTools.buildAttributedTerminiTree(splice_junctions, false, true);

		tSpliceSupportEvaluater sse = new tSpliceSupportEvaluater(strandedness, w);

		for(AnnotatedRegion r5p : junctions5p){
			Double p5p_w = (Double) r5p.getAttribute("pval_w");
			Double p5p = (Double) r5p.getAttribute("pval");
			if(p5p_w==null){
				sse.reset();
				sse.evaluate(sfr, r5p.chr, r5p.start, r5p.isNegativeStrand(), true);
				// the number of expected reads aligning is readlength + window size * depth
				p5p_w = 1.0 - Beta.regularizedBeta(1./3, sse.countIntronW.n + 1.0, sse.countExonW.n);
				p5p = 1.0 - Beta.regularizedBeta(.5, sse.countIntron.n + 1.0, sse.countExon.n);
				r5p.addAttribute("pval5p_w", p5p_w);
				r5p.addAttribute("pval5p", p5p);
				r5p.addAttribute("n5p_e_w", sse.countExonW.n);
				r5p.addAttribute("n5p_e", sse.countExon.n);
				r5p.addAttribute("n5p_i_w", sse.countIntronW.n);
				r5p.addAttribute("n5p_i", sse.countIntron.n);
			}
			if(p5p_w < alpha && p5p < 2.5*alpha){
				// if the current boundary is not involved in any junctions yet
				if(IntervalTools.BoundaryIntervals(filtered_junctions, r5p.chr, r5p.start, r5p.strand, true).size()==0){
					ExtremeObjectTracker<AnnotatedRegion, Double> j3p = new ExtremeObjectTracker<AnnotatedRegion, Double>(new Util.ComparableComparator<Double>());
					//find all the matching boundaries, tracking the most significant
					for(AnnotatedRegion jnct : IntervalTools.BoundaryIntervals(splice_junctions, r5p.chr, r5p.start, r5p.strand, true)){
						for(AnnotatedRegion r3p : junctions3p.overlappingRegions(jnct.chr, jnct.get3Prime(), jnct.get3Prime(), jnct.strand)){
							Double p3p_w = (Double) r3p.getAttribute("pval_w");
							Double p3p = (Double) r3p.getAttribute("pval");
							if(p3p==null){
								sse.reset();
								sse.evaluate(sfr, r3p.chr, r3p.start, r3p.isNegativeStrand(), false);
								p3p_w = 1.0 - Beta.regularizedBeta(1./3, sse.countIntronW.n + 1.0, sse.countExonW.n);
								p3p = 1.0 - Beta.regularizedBeta(.5, sse.countIntron.n + 1.0, sse.countExon.n);
							}
							r3p.addAttribute("pval3p_w", p3p_w);
							r3p.addAttribute("pval3p", p3p);
							r3p.addAttribute("n3p_e_w", sse.countExonW.n);
							r3p.addAttribute("n3p_e", sse.countExon.n);
							r3p.addAttribute("n3p_i_w", sse.countIntronW.n);
							r3p.addAttribute("n3p_i", sse.countIntron.n);

							j3p.put(r3p, p3p);
							if(p3p_w< alpha && p3p<2.5*alpha){
								//								System.out.println("adding "+jnct);
								//								IntervalTools.addRegion(filtered_junctions, jnct, true, false);
								if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
									filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.putAll(r5p.attributes);
									attributes.putAll(r3p.attributes);
									gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
								}
							}
						}	
					}
					// if none have significant matching boundaries, add the most significant
					if(j3p.getMin()>=alpha){
						for(AnnotatedRegion r3p : j3p.getMinObjects()){
							//							System.out.println("rescueing "+jnct);
							//							IntervalTools.addRegion(filtered_junctions, jnct, true, false);
							AnnotatedRegion jnct = IntervalTools.makeRegion(r3p.chr, r5p.start, r3p.start, r5p.isNegativeStrand());
							if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
								filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
								Map<String,Object> attributes = new HashMap<String, Object>();
								attributes.putAll(r5p.attributes);
								attributes.putAll(jnct.attributes);
								gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
						}
					}
				}			
			}
			//			else{
			//				System.out.println("rejected "+r5p);
			//			}
		}

		for(AnnotatedRegion r3p : junctions3p){
			Double p3p_w = (Double) r3p.getAttribute("pval_w");
			Double p3p = (Double) r3p.getAttribute("pval");
			if(p3p_w==null){
				sse.reset();
				sse.evaluate(sfr, r3p.chr, r3p.start, r3p.isNegativeStrand(), false);
				p3p_w = 1.0 - Beta.regularizedBeta(1./3, sse.countIntronW.n + 1.0, sse.countExonW.n);
				p3p = 1.0 - Beta.regularizedBeta(.5, sse.countIntron.n + 1.0, sse.countExon.n);
				r3p.addAttribute("pval3p_w", p3p_w);
				r3p.addAttribute("pval3p", p3p);
				r3p.addAttribute("n3p_e_w", sse.countExonW.n);
				r3p.addAttribute("n3p_e", sse.countExon.n);
				r3p.addAttribute("n3p_i_w", sse.countIntronW.n);
				r3p.addAttribute("n3p_i", sse.countIntron.n);
			}
			if(p3p_w < alpha && p3p < 2.5*alpha){
				// if the current boundary is not involved in any junctions yet
				if(IntervalTools.BoundaryIntervals(filtered_junctions, r3p.chr, r3p.start, r3p.strand, false).size()==0){
					ExtremeObjectTracker<AnnotatedRegion, Double> j5p = new ExtremeObjectTracker<AnnotatedRegion, Double>(new Util.ComparableComparator<Double>());
					//find all the matching boundaries, tracking the most significant
					for(AnnotatedRegion jnct : IntervalTools.BoundaryIntervals(splice_junctions, r3p.chr, r3p.start, r3p.strand, false)){
						for(AnnotatedRegion r5p : junctions5p.overlappingRegions(jnct.chr, jnct.get5Prime(), jnct.get5Prime(), jnct.strand)){
							Double p5p_w = (Double) r5p.getAttribute("pval_w");
							Double p5p = (Double) r5p.getAttribute("pval");
							if(p5p_w==null){
								sse.reset();
								sse.evaluate(sfr, r5p.chr, r5p.start, r5p.isNegativeStrand(), true);
								p5p_w = 1.0 - Beta.regularizedBeta(1./3, sse.countIntronW.n + 1.0, sse.countExonW.n);
								p5p = 1.0 - Beta.regularizedBeta(.5, sse.countIntron.n + 1.0, sse.countExon.n);
							}
							r5p.addAttribute("pval5p_w", p5p_w);
							r5p.addAttribute("pval5p", p5p);
							r5p.addAttribute("n5p_e_w", sse.countExonW.n);
							r5p.addAttribute("n5p_e", sse.countExon.n);
							r5p.addAttribute("n5p_i_w", sse.countIntronW.n);
							r5p.addAttribute("n5p_i", sse.countIntron.n);
							j5p.put(r5p, p5p_w);
							if(p5p_w<alpha && p5p<2.5*alpha){
								//								IntervalTools.addRegion(filtered_junctions, jnct, true, false);
								if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
									filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.putAll(r5p.attributes);
									attributes.putAll(r3p.attributes);
									gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
								}
							}
						}	
					}
					// if none have significant matching boundaries, add the most significant
					if(j5p.getMin()>=alpha){
						for(AnnotatedRegion r5p : j5p.getMinObjects()){
							//							IntervalTools.addRegion(filtered_junctions, jnct, true, false);
							AnnotatedRegion jnct = IntervalTools.makeRegion(r3p.chr, r5p.start, r3p.start, r5p.isNegativeStrand());
							if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
								filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
								Map<String,Object> attributes = new HashMap<String, Object>();
								attributes.putAll(jnct.attributes);
								attributes.putAll(r3p.attributes);
								gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
						}
					}
				}			
			}
		}
	}

	//	public static void updateJunctionCounts(SAMRecord sr, StrandedGenomicIntervalTree<Map<String,Object>> j5p, StrandedGenomicIntervalTree<Map<String,Object>> j3p){
	//		Integer lastAlignedPosition = null;
	//		int alignmentPosition = sr.getAlignmentStart();
	//		Cigar cigar = sr.getCigar();
	//
	//		for(int i=0; i<cigar.numCigarElements(); i++){
	//			CigarElement cigarElement = cigar.getCigarElement(i);
	//			if(cigarElement.getOperator().consumesReferenceBases()){
	//				boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
	//
	//				if(consumesReadBases){
	//					for(int j=0; j<cigarElement.getLength(); j++){
	//						if(lastAlignedPosition!=null){
	//							AnnotatedRegion span = new AnnotatedRegion("", sr.getReferenceName(), lastAlignedPosition, alignmentPosition, '.');
	//							if(alignmentPosition-lastAlignedPosition==1){
	//								for(char strand : Util.list('+','-')){
	//									int p5p = span.get3Prime(strand);
	//									int p3p = span.get5Prime(strand);
	//									for(AnnotatedRegion r : j5p.overlappingRegions(sr.getReferenceName(), p5p, p5p, strand)){
	//										MapCounter<String> mc = (MapCounter<String>) r.getAttribute("counts");
	//										mc.increment("n5p_i");
	//									}
	//									for(AnnotatedRegion r : j3p.overlappingRegions(sr.getReferenceName(), p3p, p3p, strand)){
	//										MapCounter<String> mc = (MapCounter<String>) r.getAttribute("counts");
	//										mc.increment("n3p_i");
	//									}
	//								}
	//								//							System.out.printf("spanned %s:%d-%d\n", sr.getReferenceName(), lastAlignedPosition, alignmentPosition);
	//							}
	//							else if(alignmentPosition-lastAlignedPosition!=1){
	//								//								System.out.printf("not spanned %s:%d-%d\n", sr.getReferenceName(), lastAlignedPosition, alignmentPosition);
	//
	//								for(char strand : Util.list('+','-')){
	//									int p5p = IntervalTools.offsetPosition(span.get5Prime(strand), 1, strand=='-', false);
	//									int p3p = IntervalTools.offsetPosition(span.get3Prime(strand), 1, strand=='-', true);
	//									for(AnnotatedRegion r : j5p.overlappingRegions(sr.getReferenceName(), p5p, p5p, strand)){
	//										MapCounter<String> mc = (MapCounter<String>) r.getAttribute("counts");
	//										mc.increment("n5p_e");
	//									}
	//									for(AnnotatedRegion r : j3p.overlappingRegions(sr.getReferenceName(), p3p, p3p, strand)){
	//										MapCounter<String> mc = (MapCounter<String>) r.getAttribute("counts");
	//										mc.increment("n3p_e");
	//									}
	//								}
	//							}
	//
	//						}
	//
	//						lastAlignedPosition = alignmentPosition;
	//						alignmentPosition++;
	//					}
	//				}
	//				else{
	//					alignmentPosition+=cigarElement.getLength();
	//				}
	//			}			
	//		}
	//	}

	public static void updateJunctionCounts(SAMRecord sr, StrandedGenomicIntervalTree<Map<String,Object>> j5p, StrandedGenomicIntervalTree<Map<String,Object>> j3p, Strandedness strandedness){
		Integer lastAlignedPosition = null;
		int alignmentPosition = sr.getAlignmentStart();
		Cigar cigar = sr.getCigar();

		StrandedGenomicIntervalSet set1_5p = new StrandedGenomicIntervalSet();
		StrandedGenomicIntervalSet set1_3p = new StrandedGenomicIntervalSet();
		StrandedGenomicIntervalSet set2_5p = new StrandedGenomicIntervalSet();
		StrandedGenomicIntervalSet set2_3p = new StrandedGenomicIntervalSet();

		int start = sr.getAlignmentStart();
		int end = sr.getAlignmentEnd();

		char[] strands = null;

		switch (strandedness) {
		case unstranded:
			strands = new char[]{'+','-'};
			break;

		default:
			strands = new char[]{BAMTools.strand(sr, strandedness)};
			break;
		}

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			if(cigarElement.getOperator().consumesReferenceBases()){
				boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();

				if(consumesReadBases){
					for(int j=0; j<cigarElement.getLength(); j++){
						if(lastAlignedPosition!=null){
							AnnotatedRegion span = new AnnotatedRegion("", sr.getReferenceName(), lastAlignedPosition, alignmentPosition, '.');
							if(alignmentPosition-lastAlignedPosition==1){

								for(char strand : strands){
									int p5p = span.get3Prime(strand);
									int p3p = span.get5Prime(strand);

									set1_5p.add(sr.getReferenceName(), p5p, p5p, strand);
									set1_3p.add(sr.getReferenceName(), p3p, p3p, strand);
								}
								//							System.out.printf("spanned %s:%d-%d\n", sr.getReferenceName(), lastAlignedPosition, alignmentPosition);
							}
							else if(alignmentPosition-lastAlignedPosition!=1){
								//								System.out.printf("not spanned %s:%d-%d\n", sr.getReferenceName(), lastAlignedPosition, alignmentPosition);

								for(char strand : strands){
									int p5p = IntervalTools.offsetPosition(span.get5Prime(strand), 1, strand=='-', false);
									int p3p = IntervalTools.offsetPosition(span.get3Prime(strand), 1, strand=='-', true);

									set2_5p.add(sr.getReferenceName(), p5p, p5p, strand);
									set2_3p.add(sr.getReferenceName(), p3p, p3p, strand);
								}
							}

						}

						lastAlignedPosition = alignmentPosition;
						alignmentPosition++;
					}
				}
				else{
					alignmentPosition+=cigarElement.getLength();
				}
			}			
		}



		for(char strand : strands){
			for(AnnotatedRegion r : j5p.overlappingRegions(sr.getReferenceName(), start, end, strand)){
				MapCounter<String> mc = (MapCounter<String>) r.getAttribute("counts");

				if(IntervalTools.isContained(set1_5p, r.chr, r.start, r.end, r.strand)){	
					mc.increment("n5p_i");	
				}
				if(IntervalTools.isContained(set2_5p, r.chr, r.start, r.end, r.strand)){	
					mc.increment("n5p_e");	
				}
			}
			for(AnnotatedRegion r : j3p.overlappingRegions(sr.getReferenceName(), start, end, strand)){
				MapCounter<String> mc = (MapCounter<String>) r.getAttribute("counts");

				if(IntervalTools.isContained(set1_3p, r.chr, r.start, r.end, r.strand)){	
					mc.increment("n3p_i");	
				}
				if(IntervalTools.isContained(set2_3p, r.chr, r.start, r.end, r.strand)){	
					mc.increment("n3p_e");	
				}
			}
		}

	}

	public static void updateJunctionCounts2(SAMRecord sr, MapFactory<Character, MapCounter<Integer>> strand_count_5p_span, MapFactory<Character, MapCounter<Integer>> strand_count_5p_splice, MapFactory<Character, MapCounter<Integer>> strand_count_3p_span, MapFactory<Character, MapCounter<Integer>> strand_count_3p_splice, Strandedness strandedness){
		Integer lastAlignedPosition = null;
		int alignmentPosition = sr.getAlignmentStart();
		Cigar cigar = sr.getCigar();

		char[] strands = null;
		// will be undefined if read is unspliced
		Character splice_strand = sr.getCharacterAttribute("XS");
		// will not be properly defined if the read is unspliced
		Boolean isNegativeStrand = splice_strand == null ? null : '-'==splice_strand;
		
		switch (strandedness) {
		case unstranded:
			strands = new char[]{'+','-'};
			break;

		default:
			strands = new char[]{BAMTools.strand(sr, strandedness)};
			break;
		}
		
		CigarElement prevCigarElement = null;
		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			
			if(cigarElement.getOperator().consumesReferenceBases()){
				boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();

				if(cigarElement.getOperator()==CigarOperator.N){
					alignmentPosition+=cigarElement.getLength();
				}
				else{
					for(int j=0; j<cigarElement.getLength(); j++){
						if(lastAlignedPosition!=null){

							if(alignmentPosition-lastAlignedPosition==1){

								for(char strand : strands){
									int p3p = IntervalTools.get3prime(lastAlignedPosition, alignmentPosition, strand=='-');
									strand_count_5p_span.get(strand).increment(p3p);
									int p5p = IntervalTools.get5prime(lastAlignedPosition, alignmentPosition, strand=='-');
									strand_count_3p_span.get(strand).increment(p5p);
								}
								//							System.out.printf("spanned %s:%d-%d\n", sr.getReferenceName(), lastAlignedPosition, alignmentPosition);
							}
							else if(alignmentPosition-lastAlignedPosition!=1 && prevCigarElement.getOperator().equals(CigarOperator.SKIPPED_REGION)){
								//								System.out.printf("not spanned %s:%d-%d\n", sr.getReferenceName(), lastAlignedPosition, alignmentPosition);

								
									int p5p = IntervalTools.offsetPosition(IntervalTools.get5prime(lastAlignedPosition, alignmentPosition, isNegativeStrand), 1, isNegativeStrand, false);
									strand_count_5p_splice.get(splice_strand).increment(p5p);
									int p3p = IntervalTools.offsetPosition(IntervalTools.get3prime(lastAlignedPosition, alignmentPosition, isNegativeStrand), 1, isNegativeStrand, true);
									strand_count_3p_splice.get(splice_strand).increment(p3p);
								
							}

						}

						lastAlignedPosition = alignmentPosition;
						alignmentPosition++;
					}
				}
			}			
			prevCigarElement = cigarElement;
		}
	}

	//	public static void countJunctionSupportingReads(SAMFileReader sfr, File splice_junction_bed, GTFWriter gw) throws FileNotFoundException{	
	//		StrandedGenomicIntervalTree<Map<String,Object>> jncts = IntervalTools.buildRegionsTree(new BEDIterator(splice_junction_bed), true, false);
	//		StrandedGenomicIntervalTree<Map<String,Object>> j5p = IntervalTools.buildAttributedTerminiTree(jncts, true, true);
	//		StrandedGenomicIntervalTree<Map<String,Object>> j3p = IntervalTools.buildAttributedTerminiTree(jncts, false, true);
	//
	//		for(AnnotatedRegion r : j5p){
	//			r.addAttribute("counts", new MapCounter<String>());
	//		}
	//
	//		for(AnnotatedRegion r : j3p){
	//			r.addAttribute("counts", new MapCounter<String>());
	//		}
	//
	//		SAMRecordIterator sri = sfr.iterator();
	//		while(sri.hasNext()){
	//			SAMRecord sr = sri.next();
	//			updateJunctionCounts(sr,j5p,j3p);
	//		}
	//		sri.close();
	//
	//		for(AnnotatedRegion r : j5p){
	//			Map<String,Object> attributes = new HashMap<String, Object>();
	//			for(String key : Util.list("n5p_e","n5p_i")){
	//				attributes.put(key, ((MapCounter<String>) r.getAttribute("counts")).get(key));
	//			}
	//			gw.write("j5p",r.chr, r.start, r.end, r.strand,AnnotatedRegion.GTFAttributeString(attributes));
	//			//			System.out.printf("%s\t%s\n", r,((MapCounter<String>) r.getAttribute("counts")).getMap());
	//		}
	//		for(AnnotatedRegion r : j3p){
	//			Map<String,Object> attributes = new HashMap<String, Object>();
	//			for(String key : Util.list("n3p_e","n3p_i")){
	//				attributes.put(key, ((MapCounter<String>) r.getAttribute("counts")).get(key));
	//			}
	//			gw.write("j3p",r.chr, r.start, r.end, r.strand,AnnotatedRegion.GTFAttributeString(attributes));
	//			//			System.out.printf("%s\t%s\n", r,((MapCounter<String>) r.getAttribute("counts")).getMap());
	//		}
	//
	//		gw.close();
	//	}

	

	public static void countJunctionSupportingReads(SAMFileReader sfr, Strandedness strandedness, GTFWriter gw) throws FileNotFoundException{	

		Instantiator<MapCounter<Integer>> i = new Instantiator<MapCounter<Integer>>() {
			public MapCounter<Integer> instantiate(Object... objects) {
				return new MapCounter<Integer>(new TreeMap<Integer,Integer>());
			}
		};

		MapFactory<Character, MapCounter<Integer>> strand_count_5p_splice = new MapFactory(i);
		MapFactory<Character, MapCounter<Integer>> strand_count_5p_span = new MapFactory(i);
		MapFactory<Character, MapCounter<Integer>> strand_count_3p_splice = new MapFactory(i);
		MapFactory<Character, MapCounter<Integer>> strand_count_3p_span = new MapFactory(i);

		SAMRecordIterator sri = sfr.iterator();
		String prev_chr = null;
		
		while(sri.hasNext()){
			
			SAMRecord sr = sri.next();
			String chr = sr.getReferenceName();
			int alignment_start = sr.getAlignmentStart();
			
			if(strandedness!=Strandedness.unstranded && sr.getAttribute("XS")!=null && BAMTools.strand(sr, strandedness)!=sr.getCharacterAttribute("XS"))
				continue;
			
			if(prev_chr!=null && !chr.equals(prev_chr)){
				writeAll(strand_count_5p_splice,strand_count_5p_span, prev_chr, true,gw);
				// IMPORTANT: corrected error when last read aligned to the chromosome is spliced
				writeAll(strand_count_3p_splice,strand_count_3p_span, prev_chr, false,gw);
				strand_count_5p_splice.getMap().clear();
				strand_count_5p_span.getMap().clear();
				strand_count_3p_splice.getMap().clear();
				strand_count_3p_span.getMap().clear();
			}
			
			//			updateJunctionCounts(sr, j5p, j3p, strandedness);
			updateJunctionCounts2(sr, strand_count_5p_span, strand_count_5p_splice, strand_count_3p_span, strand_count_3p_splice, strandedness);
			//			System.out.println(count_5p_splice);

			writeToAlignmentStart(strand_count_5p_splice, strand_count_5p_span, chr, alignment_start, true,gw);
			writeToAlignmentStart(strand_count_3p_splice, strand_count_3p_span, chr, alignment_start, false,gw);

			removeToAlignmentStart(strand_count_5p_splice,alignment_start);
			removeToAlignmentStart(strand_count_5p_span,alignment_start);
			removeToAlignmentStart(strand_count_3p_splice,alignment_start);
			removeToAlignmentStart(strand_count_3p_span,alignment_start);
			
			// IMPORTANT: corrected error when last read aligned to the chromosome is spliced
			prev_chr = chr;
		}
		sri.close();

		writeAll(strand_count_5p_splice,strand_count_5p_span, prev_chr, true,gw);
		writeAll(strand_count_3p_splice,strand_count_3p_span, prev_chr, false,gw);
	}

	public static void writeToAlignmentStart(MapFactory<Character, MapCounter<Integer>> strand_count_splice, MapFactory<Character, MapCounter<Integer>> strand_count_span, String chr, int alignment_start, boolean is5p, GTFWriter gw){
		for(char strand : new char[]{'+','-'}){
			SortedMultiIterator<Integer> smi = new SortedMultiIterator<Integer>(Integer.class, ((TreeMap<Integer,Integer>) strand_count_splice.get(strand).getMap()).navigableKeySet().iterator(), ((TreeMap<Integer,Integer>) strand_count_span.get(strand).getMap()).navigableKeySet().iterator());
			while(smi.hasNext()){
				Integer[] keys = smi.next();
				Integer key = keys[0]==null?keys[1]:keys[0];
				if(key < alignment_start){
					if(keys[0]!=null){
//						System.out.printf("writing pos %s:%d:%c %s %d %s %d\n", chr, key, strand, is5p?"5p_e":"3p_e", strand_count_3p_splice.get(strand).get(key), is5p?"5p_i":"3p_i", strand_count_3p_span.get(strand).get(key));
						String type = is5p ? "j5p" : "j3p";
						Map<String,Object> attributes = new HashMap<String, Object>();
						if(is5p){
							attributes.put("n5p_e", strand_count_splice.get(strand).get(key));
							attributes.put("n5p_i", strand_count_span.get(strand).get(key));
							}
						else{
							attributes.put("n3p_e", strand_count_splice.get(strand).get(key));
							attributes.put("n3p_i", strand_count_span.get(strand).get(key));
						}
						gw.write(type,chr, key, key, strand,AnnotatedRegion.GTFAttributeString(attributes));
					}
				}
				else{
					break;
				}
			}
		}
	}

	public static void removeToAlignmentStart(MapFactory<Character, MapCounter<Integer>> strand_count_5p_splice, int start){
		for(char strand : new char[]{'+','-'}){
			Iterator<Integer> ascendingKeySet = ((TreeMap<Integer,Integer>) strand_count_5p_splice.get(strand).getMap()).navigableKeySet().iterator();
			while(ascendingKeySet.hasNext()){
				Integer key = ascendingKeySet.next();
				//					System.out.printf("span key %d alignment start %d\n", key, sr.getAlignmentStart());
				if(key < start){
					//						System.out.printf("removing 5p_span %d %d\n", key, tree_5p_span.get(key));
					ascendingKeySet.remove();
				}
				else{
					break;
				}
			}
		}
	}

	public static void writeAll(MapFactory<Character, MapCounter<Integer>> strand_count_splice, MapFactory<Character, MapCounter<Integer>> strand_count_span, String chr, boolean is5p, GTFWriter gw){
		for(char strand : new char[]{'+','-'}){
			SortedMultiIterator<Integer> smi = new SortedMultiIterator<Integer>(Integer.class, ((TreeMap<Integer,Integer>) strand_count_splice.get(strand).getMap()).navigableKeySet().iterator(), ((TreeMap<Integer,Integer>) strand_count_span.get(strand).getMap()).navigableKeySet().iterator());
			while(smi.hasNext()){
				Integer[] keys = smi.next();
				Integer key = keys[0]==null?keys[1]:keys[0];
				if(keys[0]!=null){
//					System.out.printf("writing pos %s:%d:%c %s %d %s %d\n", chr, key, strand, is5p?"5p_e":"3p_e", strand_count_5p_splice.get(strand).get(key), is5p?"5p_i":"3p_i", strand_count_5p_span.get(strand).get(key));
					String type = is5p ? "j5p" : "j3p";
					Map<String,Object> attributes = new HashMap<String, Object>();
					if(is5p){
						attributes.put("n5p_e", strand_count_splice.get(strand).get(key));
						attributes.put("n5p_i", strand_count_span.get(strand).get(key));
						}
					else{
						attributes.put("n3p_e", strand_count_splice.get(strand).get(key));
						attributes.put("n3p_i", strand_count_span.get(strand).get(key));
					}
					gw.write(type,chr, key, key, strand,AnnotatedRegion.GTFAttributeString(attributes));
				}
			}
		}
	}

	/**
	 * 
	 * @param splice_junction_bed - bed file with all the juctions detected
	 * @param splice_count_gtf - gtf with the counts for each donor/acceptor site
	 * @param alpha - significance to accept splice junction
	 * @param gw - a list of sites that pass the filter
	 * @param gwf - the set of junctions that failed the filter
	 * @throws FileNotFoundException
	 */
	public static void filterCountedJunctions(File splice_junction_bed, File splice_count_gtf, double alpha, GTFWriter gw, GTFWriter gwf) throws FileNotFoundException{
		BEDIterator sj = new BEDIterator(splice_junction_bed);

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = IntervalTools.buildRegionsTree(sj, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> junction_counts = IntervalTools.buildRegionsTree(new TranscriptIterator(splice_count_gtf), false, true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> filtered_junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> to_filter_junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> junctions5p = IntervalTools.buildRegionsTree(new IntervalTools.AttributeIterable(junction_counts, "annotation", "j5p"), false, true);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions3p = IntervalTools.buildRegionsTree(new IntervalTools.AttributeIterable(junction_counts, "annotation", "j3p"), false, true);

		for(AnnotatedRegion r5p : junctions5p){

			int n5p_e = Integer.parseInt((String) r5p.getAttribute("n5p_e"));
			int n5p_i = Integer.parseInt((String) r5p.getAttribute("n5p_i"));

			double p5p = 1.0 - Beta.regularizedBeta(.5, n5p_i + 1.0, n5p_e);

			if(p5p < alpha){
				// if the current boundary is not involved in any junctions yet
				if(IntervalTools.BoundaryIntervals(filtered_junctions, r5p.chr, r5p.start, r5p.strand, true).size()==0){
					ExtremeObjectTracker<AnnotatedRegion, Double> j3p = new ExtremeObjectTracker<AnnotatedRegion, Double>(new Util.ComparableComparator<Double>());
					//find all the matching boundaries, tracking the most significant
					for(AnnotatedRegion jnct : IntervalTools.BoundaryIntervals(splice_junctions, r5p.chr, r5p.start, r5p.strand, true)){
						for(AnnotatedRegion r3p : junctions3p.overlappingRegions(jnct.chr, jnct.get3Prime(), jnct.get3Prime(), jnct.strand)){
							int n3p_e = Integer.parseInt((String) r3p.getAttribute("n3p_e"));
							int n3p_i = Integer.parseInt((String) r3p.getAttribute("n3p_i"));

							double p3p = 1.0 - Beta.regularizedBeta(.5, n3p_i + 1.0, n3p_e);

							j3p.put(r3p, p3p);
							if(p3p<alpha){
								//								System.out.println("adding "+jnct);
								//								IntervalTools.addRegion(filtered_junctions, jnct, true, false);
								if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
									filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.putAll(r5p.attributes);
									attributes.putAll(r3p.attributes);
									gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
								}
							}
							// if not significant, add it to the to-filter list
							else{
								if(!to_filter_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.putAll(r5p.attributes);
									attributes.putAll(r3p.attributes);
									to_filter_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand, attributes);
								}
							}
						}	
					}
					// if none have significant matching boundaries, add the most significant
					if(j3p.getMin()>=alpha){
						for(AnnotatedRegion r3p : j3p.getMinObjects()){
							//							System.out.println("rescueing "+jnct);
							//							IntervalTools.addRegion(filtered_junctions, jnct, true, false);
							AnnotatedRegion jnct = IntervalTools.makeRegion(r3p.chr, r5p.start, r3p.start, r5p.isNegativeStrand());
							if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
								filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
								Map<String,Object> attributes = new HashMap<String, Object>();
								attributes.putAll(r5p.attributes);
								attributes.putAll(r3p.attributes);
								gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
						}
					}
				}			
			}
			else{
				for(AnnotatedRegion jnct : IntervalTools.BoundaryIntervals(splice_junctions, r5p.chr, r5p.start, r5p.strand, true)){
					for(AnnotatedRegion r3p : junctions3p.overlappingRegions(jnct.chr, jnct.get3Prime(), jnct.get3Prime(), jnct.strand)){
						if(!to_filter_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
							Map<String,Object> attributes = new HashMap<String, Object>();
							attributes.putAll(r5p.attributes);
							attributes.putAll(r3p.attributes);
							to_filter_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand, attributes);
						}
					}
				}
			}
			//			else{
			//				System.out.println("rejected "+r5p);
			//			}
		}

		for(AnnotatedRegion r3p : junctions3p){
			int n3p_e = Integer.parseInt((String) r3p.getAttribute("n3p_e"));
			int n3p_i = Integer.parseInt((String) r3p.getAttribute("n3p_i"));

			double p3p = 1.0 - Beta.regularizedBeta(.5, n3p_i + 1.0, n3p_e);

			if(p3p < alpha){
				// if the current boundary is not involved in any junctions yet
				if(IntervalTools.BoundaryIntervals(filtered_junctions, r3p.chr, r3p.start, r3p.strand, false).size()==0){
					ExtremeObjectTracker<AnnotatedRegion, Double> j5p = new ExtremeObjectTracker<AnnotatedRegion, Double>(new Util.ComparableComparator<Double>());
					//find all the matching boundaries, tracking the most significant
					for(AnnotatedRegion jnct : IntervalTools.BoundaryIntervals(splice_junctions, r3p.chr, r3p.start, r3p.strand, false)){
						for(AnnotatedRegion r5p : junctions5p.overlappingRegions(jnct.chr, jnct.get5Prime(), jnct.get5Prime(), jnct.strand)){
							int n5p_e = Integer.parseInt((String) r5p.getAttribute("n5p_e"));
							int n5p_i = Integer.parseInt((String) r5p.getAttribute("n5p_i"));

							double p5p = 1.0 - Beta.regularizedBeta(.5, n5p_i + 1.0, n5p_e);

							j5p.put(r5p, p5p);
							if(p5p<alpha){
								//								IntervalTools.addRegion(filtered_junctions, jnct, true, false);
								if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
									filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.putAll(r5p.attributes);
									attributes.putAll(r3p.attributes);
									gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
								}
							}
							// if not significant, add it to the to-filter list
							else{
								if(!to_filter_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.putAll(r5p.attributes);
									attributes.putAll(r3p.attributes);
									to_filter_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
								}
							}
						}	
					}
					// if none have significant matching boundaries, add the most significant
					if(j5p.getMin()>=alpha){
						for(AnnotatedRegion r5p : j5p.getMinObjects()){
							//							IntervalTools.addRegion(filtered_junctions, jnct, true, false);
							AnnotatedRegion jnct = IntervalTools.makeRegion(r3p.chr, r5p.start, r3p.start, r5p.isNegativeStrand());
							if(!filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
								filtered_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand);
								Map<String,Object> attributes = new HashMap<String, Object>();
								attributes.putAll(r5p.attributes);
								attributes.putAll(r3p.attributes);
								gw.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(attributes));
							}
						}
					}
				}			
			}
			else{
				for(AnnotatedRegion jnct : IntervalTools.BoundaryIntervals(splice_junctions, r3p.chr, r3p.start, r3p.strand, false)){
					for(AnnotatedRegion r5p : junctions5p.overlappingRegions(jnct.chr, jnct.get5Prime(), jnct.get5Prime(), jnct.strand)){
						if(!to_filter_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
							Map<String,Object> attributes = new HashMap<String, Object>();
							attributes.putAll(r5p.attributes);
							attributes.putAll(r3p.attributes);
							to_filter_junctions.add(jnct.chr, jnct.start, jnct.end, jnct.strand, attributes);
						}
					}
				}
			}

		}

		// remove all junctions that are in the to_filter set that were rescued because they have one significant half, and were the most likely
		List<AnnotatedRegion> toRemove = new LinkedList<AnnotatedRegion>();
		for(AnnotatedRegion jnct : to_filter_junctions){
			if(filtered_junctions.contains(jnct.chr, jnct.start, jnct.end, jnct.strand)){
				toRemove.add(jnct);
			}
		}

		for(AnnotatedRegion jnct : toRemove){
			to_filter_junctions.remove(jnct.chr, jnct.start, jnct.end, jnct.strand);
		}

		for(AnnotatedRegion jnct : to_filter_junctions){
			gwf.write("splice_jnct", jnct.chr, jnct.start, jnct.end, jnct.strand, AnnotatedRegion.GTFAttributeString(jnct.attributes));
		}
	}

	public static void evaluateSpliceSupport(SAMFileReader sfr, Strandedness strandedness, File splice_junction_bed, int w, GTFWriter gw) throws FileNotFoundException{
		BEDIterator sj = new BEDIterator(splice_junction_bed);

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = IntervalTools.buildRegionsTree(sj, true, false);

		StrandedGenomicIntervalTree<Map<String,Object>> junctions5p = IntervalTools.buildAttributedTerminiTree(splice_junctions, true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> junctions3p = IntervalTools.buildAttributedTerminiTree(splice_junctions, false, true);

		tSpliceSupportEvaluater sse = new tSpliceSupportEvaluater(strandedness, w);

		for(AnnotatedRegion r5p : junctions5p){
			//			if(!r5p.chr.equals("17"))
			//				continue;
			Integer n5p_e = (Integer) r5p.getAttribute("n5p_e");

			if(n5p_e==null){
				sse.reset();
				sse.evaluate(sfr, r5p.chr, r5p.start, r5p.isNegativeStrand(), true);
				r5p.addAttribute("n5p_e_w", sse.countExonW.n);
				r5p.addAttribute("n5p_e", sse.countExon.n);
				r5p.addAttribute("n5p_i_w", sse.countIntronW.n);
				r5p.addAttribute("n5p_i", sse.countIntron.n);
			}

			gw.write("j5p", r5p.chr, r5p.start, r5p.end, r5p.strand, AnnotatedRegion.GTFAttributeString(r5p.attributes));
		}

		for(AnnotatedRegion r3p : junctions3p){
			//			if(!r3p.chr.equals("17"))
			//				continue;
			Integer n3p_e = (Integer) r3p.getAttribute("pval_w");
			if(n3p_e==null){
				sse.reset();
				sse.evaluate(sfr, r3p.chr, r3p.start, r3p.isNegativeStrand(), false);
				r3p.addAttribute("n3p_e_w", sse.countExonW.n);
				r3p.addAttribute("n3p_e", sse.countExon.n);
				r3p.addAttribute("n3p_i_w", sse.countIntronW.n);
				r3p.addAttribute("n3p_i", sse.countIntron.n);
			}

			gw.write("j3p", r3p.chr, r3p.start, r3p.end, r3p.strand, AnnotatedRegion.GTFAttributeString(r3p.attributes));
		}
	}

	public static StrandedGenomicIntervalTree<Map<String,Object>> tabulateReplicateSpliceJunctions(SAMFileReader[] sfrs, BEDWriter bw){
		Map<String,Integer> referenceLengths = BAMTools.referenceSequenceLengths(sfrs[0].getFileHeader());

		StrandedGenomicIntervalTree<Map<String,Object>> splice_junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(String chr : referenceLengths.keySet()){
			int start = 1;
			int end = referenceLengths.get(chr);
			boolean contained = true;

			for(SAMFileReader sfr : sfrs){
				SAMRecordIterator sri = sfr.query(chr, start, end, contained);
				//		List<AnnotatedRegion> splice_junctions = new LinkedList<AnnotatedRegion>();
				while(sri.hasNext()){
					SAMRecord sr = sri.next();
					if(sr.getCigarString().indexOf('N') == -1)
						continue;

					Cigar cigar = sr.getCigar();
					int alignmentPosition = sr.getAlignmentStart();

					for(int i=0; i<cigar.numCigarElements(); i++){
						CigarElement cigarElement = cigar.getCigarElement(i);
						if(cigarElement.getOperator().consumesReferenceBases()){
							boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
							if(isIntron){
								//							System.out.println(Util.list(alignmentPosition,alignmentPosition+cigarElement.getLength()-1));
								//							if((alignmentPosition >= start && alignmentPosition <= end) || (alignmentPosition+cigarElement.getLength()-1 >= start && alignmentPosition+cigarElement.getLength()-1 <= end)){
								AnnotatedRegion sj = new AnnotatedRegion("sj", chr, alignmentPosition,alignmentPosition+cigarElement.getLength()-1, sr.getCharacterAttribute("XS"));
								if(!splice_junctions.contains(sj.chr, sj.start, sj.end, sj.strand)){
									splice_junctions.add(sj);
									bw.write(sj.chr, sj.start, sj.end, sj.strand);
								}
								//							}
								//							System.out.println(alignmentPosition);
							}
							alignmentPosition += cigarElement.getLength();
							//						if(isIntron
							////								&& alignmentPosition >= start && alignmentPosition <= end
							//								){
							//							n++;
							//							System.out.println(alignmentPosition);
							//							System.out.println();
							//							break spliced;
							//						}
						}			
					}

				}
				sri.close();
			}
		}

		return splice_junctions;
	}

	public static <T extends Comparable<T>> List<T> merge(PeekingIterator<T>...lists){
		int l = lists.length;

		class IndexedIterator<T>{
			int i;
			PeekingIterator<T> pi;

			public IndexedIterator(int i, PeekingIterator<T> pi) {
				this.i = i;
				this.pi = pi;
			}
		}

		Comparator<IndexedIterator<T>> C = new Comparator<IndexedIterator<T>>() {
			public int compare(IndexedIterator<T> o1, IndexedIterator<T> o2) {
				return o1.pi.peek().compareTo(o2.pi.peek());
			}
		};

		PriorityQueue<IndexedIterator<T>> pq = new PriorityQueue<IndexedIterator<T>>(lists.length, C);
		for(int i=0; i < l; i++){
			PeekingIterator<T> q = lists[i];
			if(q.hasNext){
				pq.add(new IndexedIterator(i, q));
			}
		}

		while(!pq.isEmpty()){
			T[] array = (T[]) new Comparable[l];
			T next = null;
			do{
				IndexedIterator<T> q = pq.poll();
				next = q.pi.next();
				array[q.i] = next;
				if(q.pi.hasNext()){
					System.out.println("adding "+q.i+" back to queue after removing element "+next);
					System.out.println(q.pi.peek());
					System.out.println(q.pi.hasNext());
					pq.add(q);
				}
			} while(!pq.isEmpty() && pq.peek().pi.peek().compareTo(next)==0);
			System.out.println(Util.list(array));
		}

		return null;
	}

	static class PeekingIterator<T> implements Iterator<T>{

		Iterator<T> peekable;
		T peek;
		boolean hasNext;
		boolean hasPeek;

		public PeekingIterator(Iterator<T> peekable) {
			this.peekable = peekable;
			hasNext = peekable.hasNext();
			if(hasNext){
				peek = peekable.next();
			}
		}

		public boolean hasNext() {
			return hasNext;
		}

		public T next() {
			T next = peek;
			if(peekable.hasNext()){
				peek = peekable.next();
			}
			else{
				hasNext=false;
				peek=null;
			}
			return next;
		}

		public void remove() {
			throw new RuntimeException(this.getClass().getName()+" does not support remove() operations.");
		}

		public T peek() {
			if(!hasNext)
				throw new NoSuchElementException();
			return peek;
		}

	}

	static class SortedMultiIterator<T extends Comparable<T>> implements Iterator<T[]>{
		int l;
		PriorityQueue<IndexedIterator<T>> pq;
		Class<T> type;

		class IndexedIterator<T>{
			int i;
			PeekingIterator<T> pi;

			public IndexedIterator(int i, PeekingIterator<T> pi) {
				this.i = i;
				this.pi = pi;
			}
		}

		public SortedMultiIterator(Class<T> type, Iterator<T>...ilists) {
			this.type=type;
			l = ilists.length;

			PeekingIterator<T>[] lists = new PeekingIterator[ilists.length];
			for (int i = 0; i < lists.length; i++) {
				lists[i] = new PeekingIterator<T>(ilists[i]);
			}

			Comparator<IndexedIterator<T>> C = new Comparator<IndexedIterator<T>>() {
				public int compare(IndexedIterator<T> o1, IndexedIterator<T> o2) {
					return o1.pi.peek().compareTo(o2.pi.peek());
				}
			};

			pq = new PriorityQueue<IndexedIterator<T>>(lists.length, C);
			for(int i=0; i < l; i++){
				PeekingIterator<T> q = lists[i];
				if(q.hasNext){
					pq.add(new IndexedIterator(i, q));
				}
			}

		}

		public boolean hasNext() {
			return !pq.isEmpty();
		}

		public T[] next() {

			T[] array = (T[]) Array.newInstance(type, l);
			T next = null;
			do{
				IndexedIterator<T> q = pq.poll();
				next = q.pi.next();
				array[q.i] = next;
				if(q.pi.hasNext()){
					//				System.out.println("adding "+q.i+" back to queue after removing element "+next);
					//				System.out.println(q.pi.peek());
					//				System.out.println(q.pi.hasNext());
					pq.add(q);
				}
			} while(!pq.isEmpty() && pq.peek().pi.peek().compareTo(next)==0);

			return array;
		}

		public void remove() {
			throw new RuntimeException(this.getClass().getName()+" does not support remove() operation");
		}

	}

	public static void main(String[] args) throws FileNotFoundException {


		if(true){
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/SRP017778/merged/SRR645855.bam"));
			System.out.println(BAMTools.totalAlignedReads(sfr));
			//			countJunctionSupportingReads(sfr, new File("/mnt/LaiLab/sol/GSE51572/test/isoscm/tmp/g.sj.bed"), new GTFWriter("/mnt/LaiLab/sol/GSE51572/test/g1.counts.gtf"));
			countJunctionSupportingReads(sfr, Strandedness.unstranded, new GTFWriter(System.out));
		}
		if(false){
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			File splice_junction_bed = new File("/home/sol/lailab/sol/GSE41637/nb/gtf/tmp/SRR594393.sj.bed");
			//			GTFWriter gw = new GTFWriter("/dev/stdout");
			//			tFilterSpliceJunctions(sfr, Strandedness.reverse_forward, splice_junction_bed, 50, .05, gw);
			GTFWriter gw = new GTFWriter("SRR594393.sj.all.counted.gtf");
			evaluateSpliceSupport(sfr, Strandedness.reverse_forward, splice_junction_bed, 50, gw);
			gw.close();
			//			GTFWriter filtered = new GTFWriter("SRR594393.sj.all.filtered.gtf");
			//			filterCountedJunctions(splice_junction_bed, new File("SRR594393.sj.17.counted.gtf"), .05, filtered);
			//			filtered.close();

		}
		if(false)
		{
			SpliceSupportEvaluater sse = new SpliceSupportEvaluater(Strandedness.reverse_forward, 5);
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/piero/Project_3949/mapped/YH.bam"));
			sse.evaluate(sfr, "v2_chr3L_random_348", 6028, true, false);
			System.out.println(sse.countExon.n);
			System.out.println(sse.countSpanning.n);
			System.out.println();
			sse.reset();
			sse.evaluate(sfr, "v2_chr3L_random_348", 6081, true, true);
			System.out.println(sse.countExon.n);
			System.out.println(sse.countSpanning.n);
			System.out.println();
		}
		if(false)
		{
			//			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			//			File splice_junction_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.sj.bed");
			//			BEDWriter bw = new BEDWriter("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.sj.filtered.bed");
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/workspace/IsoSCM/tests/enah.brain.bam"));
			File splice_junction_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/enah.sj.bed");
			BEDWriter bw = new BEDWriter("/home/sol/workspace/IsoSCM/tests/gtf/tmp/enah.sj.filtered.bed");

			filterSpliceJunctions(sfr, Strandedness.reverse_forward, splice_junction_bed, 5, .5, bw);
			bw.close();
		}
		if(false){
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/data/sangercenter/hippocampus.bam"));

			//		spliceJunction(sfr, new BEDWriter(System.out));


			System.out.println(spliceJunction(sfr, "chr1", 196383823, 196617248, false));

			StrandedGenomicIntervalTree<Map<String,Object>> junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();
			for(AnnotatedRegion junction : spliceJunction(sfr, "chr1", 196383823, 196617248, false)){
				boolean contained = false;
				for(AnnotatedRegion overlap : junctions.overlappingRegions(junction.chr, junction.start, junction.end, junction.strand)){
					if(overlap.start == junction.start && overlap.end == junction.end){
						contained=true;
					}
				}
				if(!contained){
					System.out.println(junction);
					junctions.add(junction);

				}
			}
		}
	}
}
