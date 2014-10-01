package scm;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import multisample.JointSegmentationResult;
import net.sf.samtools.SAMFileReader;

import org.apache.commons.math3.special.Beta;

import changepoint.ChangePoint;
import tools.AnnotatedRegion;
import tools.BAMTools;
import tools.GTFTools.GTFWriter;
import tools.IntervalTools;
import tools.ParseGTF.TranscriptIterator;
import tools.SegmentFragmenter;
import tools.StrandedGenomicIntervalSet;
import tools.StrandedGenomicIntervalTree;
import tools.Strandedness;
import util.Util;
import util.Util.ExtremeTracker;

public class IdentifyChangePoints {

	public static StrandedGenomicIntervalTree<Map<String, Object>> identifyUnconstrainedNegativeBinomialPoints(SAMFileReader sfr, String chr, int start, int end, int maxBins, int binSize, int minCP, Strandedness strandedness, boolean isNegativeStrand, double alpha_0, double beta_0, int nb_r, int r, double p, double min_fold){

		int l = end-start+1;

		int nBins = (l+binSize-1)/binSize;

		// identify change points
		// cluster them, selecting the most likely from each cluster
		// write change points to GTF

		StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();

		char strand = isNegativeStrand ? '-' : '+';

		if(nBins>maxBins){
			nBins=maxBins;
			binSize = (l+maxBins-1)/maxBins;
			System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
		}


		// calculate the coverage in this chunk
		int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
		double[] y = new double[nBins];

		for(int i=0; i<nBins; i++){
			y[i] = cov[i];
		}

		//		SegmentationResult segmentation = BetaNegativeBinomial.viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p);
		SegmentationResult segmentation = IdentifyConstrainedChangePoints.fold_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, min_fold);

		StrandedGenomicIntervalTree<Map<String,Object>> gw = new StrandedGenomicIntervalTree<Map<String,Object>>();
		if(segmentation.n_segments>1){
			//				int[][] segments = segmentation.segments();
			double[] segment_mle = segmentation.segment_mle;
			int[] changepoints = segmentation.change_point;

			for(int i=0 ; i<changepoints.length; i++){

				int changepoint_i = changepoints[i];
				// if it's an increase in coverage, the annotation should be shifted over by one bin
				if(segment_mle[i]<segment_mle[i+1]){
					changepoint_i++;
				}

				if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

					int position = start + (changepoint_i*binSize) + (binSize/2);

					if(position >= start && position <= end){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("before_mle", segment_mle[i]);
						attributes.put("after_mle", segment_mle[i+1]);
						attributes.put("type","changepoint");

						changePoints.add(chr, position, position, strand, attributes);
						// write each changepoint 
						gw.add(chr, position, position, strand, attributes);
					}
				}
			}
		}

		return gw;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> identifyConstrainedNegativeBinomialPoints(SAMFileReader sfr, String chr, int start, int end, int maxBins, int binSize, int minCP, Strandedness strandedness, boolean isNegativeStrand, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing, double min_fold){


		int l = end-start+1;

		int nBins = (l+binSize-1)/binSize;

		char strand = isNegativeStrand ? '-' : '+';

		//TODO adjust this and unconstrained method for maxbins situation

		if(nBins>maxBins){
			nBins=maxBins;
			binSize = (l+maxBins-1)/maxBins;
			System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
		}

		// identify change points
		// cluster them, selecting the most likely from each cluster
		// write change points to GTF

		StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();



		// calculate the coverage in this chunk
		int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
		double[] y = new double[nBins];

		for(int i=0; i<nBins; i++){
			y[i] = cov[i];
		}

		//		SegmentationResult segmentation = TestConstraint.constrained_viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing);
		SegmentationResult segmentation = IdentifyConstrainedChangePoints.doubly_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing, min_fold);

		if(segmentation.n_segments>1){
			//				int[][] segments = segmentation.segments();
			double[] segment_mle = segmentation.segment_mle;
			int[] changepoints = segmentation.change_point;

			for(int i=0 ; i<changepoints.length; i++){

				int changepoint_i = changepoints[i];
				// if it's an increase in coverage, the annotation should be shifted over by one bin
				if(segment_mle[i]<segment_mle[i+1]){
					changepoint_i++;					
				}

				if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

					int position = start + (changepoint_i*binSize) + (binSize/2);

					if(!(segment_mle[i]>segment_mle[i+1] == constrain_decreasing))
						throw new RuntimeException(Util.sprintf("Unexpected changepoint %s:%d:%c",chr,position,strand));

					// make sure that the binning did not result in a cp outside of the queried segment
					if(position >= start && position <= end){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("before_mle", segment_mle[i]);
						attributes.put("after_mle", segment_mle[i+1]);
						attributes.put("type","changepoint");

						changePoints.add(chr, position, position, strand, attributes);
					}
				}
			}
		}

		return changePoints;
	}

	public static void identifyConstrainedNegativeBinomialPoints(SAMFileReader sfr, String chr, int start, int end, int maxBins, int binSize, int minCP, Strandedness strandedness, boolean isNegativeStrand, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing, double min_fold, GTFWriter gw){


		int l = end-start+1;

		int nBins = (l+binSize-1)/binSize;

		char strand = isNegativeStrand ? '-' : '+';

		//TODO adjust this and unconstrained method for maxbins situation

		if(nBins>maxBins){
			nBins=maxBins;
			binSize = (l+maxBins-1)/maxBins;
			System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
		}

		// identify change points
		// cluster them, selecting the most likely from each cluster
		// write change points to GTF

		StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();



		// calculate the coverage in this chunk
		int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
		double[] y = new double[nBins];

		for(int i=0; i<nBins; i++){
			y[i] = cov[i];
		}

		//		SegmentationResult segmentation = TestConstraint.constrained_viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing);
		SegmentationResult segmentation = IdentifyConstrainedChangePoints.doubly_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing, min_fold);

		if(segmentation.n_segments>1){
			//				int[][] segments = segmentation.segments();
			double[] segment_mle = segmentation.segment_mle;
			int[] changepoints = segmentation.change_point;

			for(int i=0 ; i<changepoints.length; i++){

				int changepoint_i = changepoints[i];
				// if it's an increase in coverage, the annotation should be shifted over by one bin
				if(segment_mle[i]<segment_mle[i+1]){
					changepoint_i++;					
				}

				if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

					int position = start + (changepoint_i*binSize) + (binSize/2);

					if(!(segment_mle[i]>segment_mle[i+1] == constrain_decreasing))
						throw new RuntimeException(Util.sprintf("Unexpected changepoint %s:%d:%c",chr,position,strand));

					// make sure that the binning did not result in a cp outside of the queried segment
					if(position > start && position < end){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("before_mle", segment_mle[i]);
						attributes.put("after_mle", segment_mle[i+1]);
						attributes.put("type","changepoint");

						changePoints.add(chr, position, position, strand, attributes);
						// write each changepoint 
						gw.write("exon",chr, position, position, strand, AnnotatedRegion.GTFAttributeString(attributes));
					}
				}
			}
		}
	}

	public static void identifyAConstrainedNegativeBinomialPoints(SAMFileReader sfr, String chr, int start, int end, int maxBins, int binSize, int minCP, Strandedness strandedness, boolean isNegativeStrand, double alpha_0, double beta_0, int nb_r, int r, double p, double min_fold, GTFWriter gw){


		int l = end-start+1;

		int nBins = (l+binSize-1)/binSize;

		char strand = isNegativeStrand ? '-' : '+';

		//TODO adjust this and unconstrained method for maxbins situation

		if(nBins>maxBins){
			nBins=maxBins;
			binSize = (l+maxBins-1)/maxBins;
			System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
		}

		// identify change points
		// cluster them, selecting the most likely from each cluster
		// write change points to GTF

		StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();

		// calculate the coverage in this chunk
		int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
		double[] y = new double[nBins];

		for(int i=0; i<nBins; i++){
			y[i] = cov[i];
		}

		//		SegmentationResult segmentation = TestConstraint.constrained_viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing);
		SegmentationResult segmentation = IdentifyConstrainedChangePoints.a_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, min_fold);

		if(segmentation.n_segments>1){
			//				int[][] segments = segmentation.segments();
			double[] segment_mle = segmentation.segment_mle;
			int[] changepoints = segmentation.change_point;

			for(int i=0 ; i<changepoints.length; i++){

				int changepoint_i = changepoints[i];
				// if it's an increase in coverage, the annotation should be shifted over by one bin
				if(segment_mle[i]<segment_mle[i+1]){
					changepoint_i++;					
				}

				if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

					int position = start + (changepoint_i*binSize) + (binSize/2);

					// make sure that the binning did not result in a cp outside of the queried segment
					if(position > start && position < end){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("before_mle", segment_mle[i]);
						attributes.put("after_mle", segment_mle[i+1]);
						attributes.put("type","changepoint");

						changePoints.add(chr, position, position, strand, attributes);
						// write each changepoint 
						gw.write("exon",chr, position, position, strand, AnnotatedRegion.GTFAttributeString(attributes));
					}
				}
			}
		}
	}

	public static List<ChangePoint> identifyConstrainedNegativeBinomialPoints(String[] ids, SAMFileReader[] sfrs, Strandedness[] strandednesses, String chr, int start, int end, int maxBins, int binSize, int minCP, boolean isNegativeStrand, double alpha_0, double beta_0, int nb_r, int r, double p, boolean constrain_decreasing, double min_fold){


		int l = end-start+1;

		int nBins = (l+binSize-1)/binSize;

		char strand = isNegativeStrand ? '-' : '+';

		//TODO adjust this and unconstrained method for maxbins situation

		if(nBins>maxBins){
			nBins=maxBins;
			binSize = (l+maxBins-1)/maxBins;
			System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
		}

		// identify change points
		// cluster them, selecting the most likely from each cluster
		// write change points to GTF

		List<ChangePoint> changePoints = new ArrayList<ChangePoint>();



		// calculate the coverage in this chunk
		double[][] y = new double[sfrs.length][];
		for(int i=0; i<sfrs.length; i++){
			int[] cov = BAMTools.binnedMaxEndCoverage(sfrs[i], chr, start, nBins, binSize, false, true, strandednesses[i], isNegativeStrand);
			y[i] = new double[nBins];
			for(int j=0; j<nBins; j++){
				y[i][j] = cov[j];
			}
		}

		//		SegmentationResult segmentation = TestConstraint.constrained_viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing);
		JointSegmentationResult segmentation = IdentifyConstrainedChangePoints.doubly_constrained_multi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing, min_fold);

		if(segmentation.n_segments>1){
			//				int[][] segments = segmentation.segments();
			double[][] segment_mle = segmentation.segment_mle;
			int[] changepoints = segmentation.change_point;

			List<Integer> positions = new ArrayList<Integer>(segmentation.n_segments);
			List<double[]> before_mles = new ArrayList<double[]>(segmentation.n_segments);
			List<double[]> after_mles = new ArrayList<double[]>(segmentation.n_segments);
			
			
			StrandedGenomicIntervalTree<Map<String,Object>> segments = new StrandedGenomicIntervalTree<Map<String,Object>>();
			segments.add(chr,start,end,strand);
			
			for(int i=0 ; i<changepoints.length; i++){

				int changepoint_i = changepoints[i];

				if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

					int position = start + (changepoint_i*binSize) + (binSize/2);

					// make sure that the binning did not result in a cp outside of the queried segment
					if(position >= start && position <= end){
						positions.add(position);
						
						before_mles.add(segment_mle[i]);
						after_mles.add(segment_mle[i+1]);

						SegmentFragmenter.fragment(segments, chr, position, strand, isNegativeStrand);
					}
				}
			}
			
			List<AnnotatedRegion> slist = new ArrayList<AnnotatedRegion>();
			for(AnnotatedRegion s : segments){
				slist.add(s);
			}
				
			for(int i=0; i<positions.size(); i++){
				AnnotatedRegion before = slist.get(i);
				AnnotatedRegion after = slist.get(i+1);
				AnnotatedRegion pos = new AnnotatedRegion("changepoint",chr,positions.get(i),positions.get(i),strand);
				double[] before_mle = before_mles.get(i);
				double[] after_mle = after_mles.get(i);
				
				changePoints.add(new ChangePoint(pos, isNegativeStrand?after:before, !isNegativeStrand?after:before, ids, isNegativeStrand?after_mle:before_mle, !isNegativeStrand?after_mle:before_mle));
			}
			
			
		}

		return changePoints;
	}

	public static void identifyUnconstrainedNegativeBinomialPoints(SAMFileReader sfr, String chr, int start, int end, int maxBins, int binSize, int minCP, Strandedness strandedness, boolean isNegativeStrand, double alpha_0, double beta_0, int nb_r, int r, double p, double min_fold, GTFWriter gw){

		int l = end-start+1;

		int nBins = (l+binSize-1)/binSize;

		// identify change points
		// cluster them, selecting the most likely from each cluster
		// write change points to GTF

		StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();

		char strand = isNegativeStrand ? '-' : '+';

		if(nBins>maxBins){
			nBins=maxBins;
			binSize = (l+maxBins-1)/maxBins;
			System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
		}


		// calculate the coverage in this chunk
		int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
		double[] y = new double[nBins];

		for(int i=0; i<nBins; i++){
			y[i] = cov[i];
		}

		//		SegmentationResult segmentation = BetaNegativeBinomial.viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p);
		SegmentationResult segmentation = IdentifyConstrainedChangePoints.fold_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, min_fold);

		if(segmentation.n_segments>1){
			//				int[][] segments = segmentation.segments();
			double[] segment_mle = segmentation.segment_mle;
			int[] changepoints = segmentation.change_point;

			for(int i=0 ; i<changepoints.length; i++){

				int changepoint_i = changepoints[i];
				// if it's an increase in coverage, the annotation should be shifted over by one bin
				if(segment_mle[i]<segment_mle[i+1]){
					changepoint_i++;
				}

				if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

					int position = start + (changepoint_i*binSize) + (binSize/2);

					if(position > start && position < end){
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("before_mle", segment_mle[i]);
						attributes.put("after_mle", segment_mle[i+1]);
						attributes.put("type","changepoint");

						changePoints.add(chr, position, position, strand, attributes);
						// write each changepoint 
						gw.write("exon",chr, position, position, strand, AnnotatedRegion.GTFAttributeString(attributes));
					}
				}
			}
		}

	}

	public static void calculate_spliced_coverage(SAMFileReader sfr, File assembly_gtf, Strandedness strandedness, GTFWriter coverage_gtf_writer) throws FileNotFoundException{
		StrandedGenomicIntervalTree<Map<String, Object>> continuousSegments = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();

		// could do one chromosome at a time to avoid loading all annotations into memory
		TranscriptIterator ti = new TranscriptIterator(assembly_gtf);
		for(AnnotatedRegion segment : ti){
			if("exon".equals(segment.annotation)){
				SegmentFragmenter.add(continuousSegments, segment.chr, segment.start, segment.end, segment.strand);		
				exons.add(segment.chr, segment.start, segment.end, segment.strand, segment.attributes);
			}
		}

		StrandedGenomicIntervalTree<Map<String, Object>> quantified = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion segment : continuousSegments){
			Double med = null;
			if(strandedness==Strandedness.unstranded){
				med = Util.median(Util.list(BAMTools.coverage(sfr, segment.chr, segment.start, segment.end, false, strandedness)[0]));
			}
			else{
				med = Util.median(Util.list(BAMTools.coverage(sfr, segment.chr, segment.start, segment.end, false, strandedness)[segment.isNegativeStrand()?1:0]));
			}
			//			int nc = BAMTools.nConsumingReads(sfr, segment.chr, segment.start, segment.end, false, strandedness, segment.isNegativeStrand());
			//			System.out.println(segment.attributes);
			Map<String,Object> attributes = new HashMap<String, Object>();
			attributes.put("coverage", med);
			quantified.add(segment.chr, segment.start, segment.end, segment.strand, attributes);
			//			if(quantified.size()>100)
			//				break;
			//			if(med>10)
			//				System.out.printf("%s\t%.1f\t%d\n", segment, med, nc);
			// add attributes to the segments
		}

		ti = new TranscriptIterator(assembly_gtf);
		for(AnnotatedRegion e : ti){
			if("exon".equals(e.annotation)){		
				//			List<Double> l = new ArrayList<Double>();
				//			List<AnnotatedRegion> ll = new ArrayList<AnnotatedRegion>();
				ExtremeTracker<Double> min_coverage=new ExtremeTracker<Double>();
				for(AnnotatedRegion r : IntervalTools.ContainedIntervals(quantified, e.chr, e.start, e.end, e.strand)){
					//				l.add((Double) r.getAttribute("coverage"));
					//				ll.add(r);
					min_coverage.put((Double) r.getAttribute("coverage"));
				}
				//			chr19:3936269-3938312
				//			if(l.size()>0)
				//			System.out.printf("%s\t%s\t%s\n", e,l,ll);

				if(min_coverage.hasExtrema){
					e.addAttribute("coverage", min_coverage.getMin());		
				}

				//			}
			}

			coverage_gtf_writer.write(e.annotation, e.chr, e.start, e.end, e.strand, AnnotatedRegion.GTFAttributeString(e.attributes));
		}

		// go through the exons, select the overlapping segments, attach the minimum coverage of any overlapping segment.
	}

	public static void calculate_coverage(SAMFileReader sfr, File unspliced_gtf, Strandedness strandedness, GTFWriter coverage_gtf_writer) throws FileNotFoundException{
		//		StrandedGenomicIntervalTree<Map<String, Object>> continuousSegments = IntervalTools.buildRegionsTree(new TranscriptIterator(unspliced_gtf), true, true);
		TranscriptIterator ti =new TranscriptIterator(unspliced_gtf);

		for(AnnotatedRegion segment : ti){
			//			System.out.println(segment.chr);
			//			if(e.chr.equals("3")){
			//			System.out.println(segment);
			//			try{
			//				System.out.println(BAMTools.nAlignedReads(sfr, segment, false));

			Double med = null;
			if(strandedness==Strandedness.unstranded){
				med = Util.median(Util.list(BAMTools.coverage(sfr, segment.chr, segment.start, segment.end, false, strandedness)[0]));
			}
			else{
				med = Util.median(Util.list(BAMTools.coverage(sfr, segment.chr, segment.start, segment.end, false, strandedness)[segment.isNegativeStrand()?1:0]));
			}
			segment.addAttribute("coverage", med);
			//				int nc = BAMTools.nConsumingReads(sfr, segment.chr, segment.start, segment.end, false, strandedness, segment.isNegativeStrand());
			//				if(med>10)
			//				System.out.printf("%s\t%.1f\n", e, med);
			//				System.out.printf("%s\t%.1f\t%d\n", segment, med, nc);
			//			}
			//			catch(Exception e){
			//				e.printStackTrace();
			//			}
			coverage_gtf_writer.write(segment.annotation, segment.chr, segment.start, segment.end, segment.strand, AnnotatedRegion.GTFAttributeString(segment.attributes));	
			//			}
		}
	}
	
	

	/* TODO
	public static void calculate_coverage(BBFileReader pos_cov, File unspliced_gtf, Strandedness strandedness, GTFWriter coverage_gtf_writer) throws FileNotFoundException{
		//		StrandedGenomicIntervalTree<Map<String, Object>> continuousSegments = IntervalTools.buildRegionsTree(new TranscriptIterator(unspliced_gtf), true, true);
		TranscriptIterator ti =new TranscriptIterator(unspliced_gtf);

		for(AnnotatedRegion segment : ti){
			//			System.out.println(segment.chr);
			//			if(e.chr.equals("3")){
			//			System.out.println(segment);
			//			try{
			//				System.out.println(BAMTools.nAlignedReads(sfr, segment, false));

			Double med = null;
			if(strandedness==Strandedness.unstranded){
				med = Util.median(Util.list(BAMTools.coverage(sfr, segment.chr, segment.start, segment.end, false, strandedness)[0]));
			}
			else{
				med = Util.median(Util.list(BAMTools.coverage(sfr, segment.chr, segment.start, segment.end, false, strandedness)[segment.isNegativeStrand()?1:0]));
			}
			segment.addAttribute("coverage", med);
			//				int nc = BAMTools.nConsumingReads(sfr, segment.chr, segment.start, segment.end, false, strandedness, segment.isNegativeStrand());
			//				if(med>10)
			//				System.out.printf("%s\t%.1f\n", e, med);
			//				System.out.printf("%s\t%.1f\t%d\n", segment, med, nc);
			//			}
			//			catch(Exception e){
			//				e.printStackTrace();
			//			}
			coverage_gtf_writer.write(segment.annotation, segment.chr, segment.start, segment.end, segment.strand, AnnotatedRegion.GTFAttributeString(segment.attributes));	
			//			}
		}
	}
	*/

	public static void identifyNegativeBinomialChangePointsInLongSegments(SAMFileReader sfr, File spliced_exon_gtf, File intronic_exon_gtf, GTFWriter changepointWriter, int minLength, int maxBins, int binSize, int minCP, Strandedness strandedness, double alpha_0, double beta_0, int nb_r, int r, double p, boolean internal, double min_fold, int min_terminal) throws FileNotFoundException{
		StrandedGenomicIntervalTree<Map<String, Object>> continuousSegments = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(File gtf : Util.list(spliced_exon_gtf, intronic_exon_gtf)){
			TranscriptIterator ti = new TranscriptIterator(gtf);
			for(AnnotatedRegion segment : ti){
				SegmentFragmenter.add(continuousSegments, segment.chr, segment.start, segment.end, segment.strand);		
				exons.add(segment.chr, segment.start, segment.end, segment.strand, segment.attributes);
			}
		}

		for(AnnotatedRegion segment : continuousSegments){
			//			System.out.printf("segment %s\n", segment);
			if(segment.size-min_terminal >= minLength){

				Set<String> annotations = new HashSet<String>();
				for(AnnotatedRegion exon : exons.overlappingRegions(segment.chr, segment.start, segment.end, segment.strand)){
					String type = (String) exon.getAttribute("type");
					if(type!=null){
						annotations.add(type);
					}
				}

				if(annotations.size()==1){
					String type = annotations.iterator().next();
					if(internal && "internal_exon".equals(type)){
						identifyUnconstrainedNegativeBinomialPoints(sfr, segment.chr, segment.start, segment.end, maxBins, binSize, minCP, strandedness, segment.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, min_fold, changepointWriter);
					}
					else if("5p_exon".equals(type)){
						boolean constrained_decreasing = segment.isNegativeStrand();
						int[] i = IntervalTools.offsetInterval(segment.start, segment.end, -min_terminal, 0, segment.isNegativeStrand());
						identifyConstrainedNegativeBinomialPoints(sfr, segment.chr, i[0], i[1], maxBins, binSize, minCP, strandedness, segment.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, constrained_decreasing, min_fold, changepointWriter);
					}
					else if("3p_exon".equals(type)){
						boolean constrained_decreasing = !segment.isNegativeStrand();
						int[] i = IntervalTools.offsetInterval(segment.start, segment.end, 0, -min_terminal, segment.isNegativeStrand());
						identifyConstrainedNegativeBinomialPoints(sfr, segment.chr, i[0], i[1], maxBins, binSize, minCP, strandedness, segment.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, constrained_decreasing, min_fold, changepointWriter);
					}
				}
			}
		}
	}

	public static void identifyNegativeBinomialChangePointsInUnsplicedSegments(SAMFileReader sfr, File unspliced_gtf, File spliced_gtf, GTFWriter changepointWriter, int minLength, int maxBins, int binSize, int minCP, Strandedness strandedness, double alpha_0, double beta_0, int nb_r, int r, double p, double min_fold, int min_terminal) throws FileNotFoundException{

		StrandedGenomicIntervalTree<Map<String,Object>> spliced_exon_boundaries = IntervalTools.buildTerminiTree(new TranscriptIterator(spliced_gtf),true,false);
		TranscriptIterator ti = new TranscriptIterator(unspliced_gtf);

		for(AnnotatedRegion segment : ti){
			//			System.out.printf("segment %s\n", segment);
			if(segment.size-min_terminal >= minLength){

				if(!spliced_exon_boundaries.overlappingRegions(segment.chr, segment.start-1, segment.end+1, segment.strand).iterator().hasNext()){	
					int[] i = IntervalTools.offsetInterval(segment.start, segment.end, -min_terminal, -min_terminal, segment.isNegativeStrand());
					identifyAConstrainedNegativeBinomialPoints(sfr, segment.chr, i[0], i[1], maxBins, binSize, minCP, strandedness, segment.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, min_fold, changepointWriter);					
				}
			}
		}
	}

	public static void identifyComposite2() throws FileNotFoundException{


		StrandedGenomicIntervalTree<Map<String,Object>> junction_counts = IntervalTools.buildRegionsTree(new TranscriptIterator("SRR594393.sj.all.counted.gtf"), false, true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> junctions5p = IntervalTools.buildRegionsTree(new IntervalTools.AttributeIterable(junction_counts, "annotation", "j5p"), false, true);

		TranscriptIterator tiu = new TranscriptIterator("/home/sol/lailab/sol/GSE41637/300n/SRR594393.unspliced.gtf");
		TranscriptIterator ti = new TranscriptIterator("/home/sol/lailab/sol/GSE41637/300n/SRR594393.coverage.gtf");
		
		GTFWriter gw = new GTFWriter("SRR594393.composite.gtf");
		
		StrandedGenomicIntervalSet unspliced = IntervalTools.buildStrandedIntervalSet(tiu);
		StrandedGenomicIntervalTree<Map<String,Object>> isoscm = IntervalTools.buildRegionsTree(ti, false, true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> internal = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> spliced_exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> splice_jnct = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> splice_segs = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : isoscm){
			if("exon".equals(r.getAttribute("annotation"))){
				spliced_exons.add(r.chr, r.start, r.end, r.strand, r.attributes);
				if("internal_exon".equals(r.getAttribute("type"))){
					internal.add(r.chr, r.start, r.end, r.strand, r.attributes);
				}
			}
			if("splice_jnct".equals(r.getAttribute("annotation"))){
				splice_jnct.add(r.chr, r.start, r.end, r.strand, r.attributes);
				SegmentFragmenter.add(splice_segs, r.chr, r.start, r.end, r.strand);
			}
		}

		StrandedGenomicIntervalTree<Map<String,Object>> internal3p = IntervalTools.buildTerminiTree(internal, false, true, false);

		SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
		double alpha = .05;

		for(AnnotatedRegion e : internal3p){
			List<AnnotatedRegion> splices = IntervalTools.BoundaryIntervals(splice_jnct, e.chr, IntervalTools.offsetPosition(e.get3Prime(),1,e.isNegativeStrand(),false), e.strand, true);
			List<AnnotatedRegion> segs = IntervalTools.BoundaryIntervals(splice_segs, e.chr, IntervalTools.offsetPosition(e.get3Prime(),1,e.isNegativeStrand(),false), e.strand, true);
			if(splices.size()>0 && segs.size()>0){
				Collections.sort(splices,new Comparator<AnnotatedRegion>() {
					public int compare(AnnotatedRegion o1, AnnotatedRegion o2) {
						return o1.size-o2.size;
					}
				});

				AnnotatedRegion shortest = splices.get(0);
				AnnotatedRegion seg = segs.get(0);
				List<AnnotatedRegion> un = IntervalTools.OverlappingIntervals(unspliced, e.chr, seg.start, seg.end, seg.strand);
				if(seg.size > 500 && IntervalTools.OverlappingIntervals(spliced_exons, e.chr, seg.start, seg.end, seg.strand).size()==0 && un.size()==1){
					AnnotatedRegion us = un.get(0); 

					int pos = IntervalTools.offsetPosition(e.get3Prime(), 1, e.isNegativeStrand(), false);

					Iterator<AnnotatedRegion> it = junctions5p.overlappingRegions(e.chr, pos, pos, e.strand).iterator();
					if(it.hasNext()){
						AnnotatedRegion r5p = it.next();

						int n5p_e_w = Integer.parseInt((String) r5p.getAttribute("n5p_e_w"));
						int n5p_e = Integer.parseInt((String) r5p.getAttribute("n5p_e"));
						int n5p_i_w = Integer.parseInt((String) r5p.getAttribute("n5p_i_w"));
						int n5p_i = Integer.parseInt((String) r5p.getAttribute("n5p_i"));

						double p5p_w = 1.0 - Beta.regularizedBeta(1./3, n5p_i_w + 1.0, n5p_e_w);
						double p5p = 1.0 - Beta.regularizedBeta(.5, n5p_i + 1.0, n5p_e);
						double pBp = Beta.regularizedBeta(.3/3, n5p_i_w + 1.0, n5p_e_w);

						if(p5p_w < alpha && p5p < 2.5*alpha && pBp < alpha){
//							System.out.printf("seg %s\n", seg);
							
							String chr = seg.chr;
							int start = seg.start;
							int end = seg.end;
							int binSize=20;
							char strand = seg.strand;
							boolean isNegativeStrand=seg.isNegativeStrand();
							int maxBins=2000;
							Strandedness strandedness = Strandedness.reverse_forward;
							double alpha_0=1; double beta_0=1;
							int minCP=1;
							double nb_r=5;
							int r=10;
							double p=.99;
							//			boolean constrain_decreasing=false;
							double min_fold=.333;

							int l = end-start+1;

							int nBins = (l+binSize-1)/binSize;

							//			char strand = isNegativeStrand ? '-' : '+';


							if(nBins>maxBins){
								nBins=maxBins;
								binSize = (l+maxBins-1)/maxBins;
								System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
							}

							// identify change points
							// cluster them, selecting the most likely from each cluster
							// write change points to GTF

							StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();



							// calculate the coverage in this chunk
							int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
							double[] y = new double[nBins];

							for(int i=0; i<nBins; i++){
								y[i] = cov[i];
							}

							//		SegmentationResult segmentation = TestConstraint.constrained_viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing);
							SegmentationResult segmentation = IdentifyConstrainedChangePoints.doubly_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, !seg.isNegativeStrand(), min_fold);

							if(segmentation.n_segments>1){
								//				int[][] segments = segmentation.segments();
								double[] segment_mle = segmentation.segment_mle;
								int[] changepoints = segmentation.change_point;

								for(int i=0 ; i<changepoints.length; i++){

									int changepoint_i = changepoints[i];
									// if it's an increase in coverage, the annotation should be shifted over by one bin
									if(segment_mle[i]<segment_mle[i+1]){
										changepoint_i++;					
									}

									if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

										int position = start + (changepoint_i*binSize) + (binSize/2);

										//						if(!(segment_mle[i]>segment_mle[i+1] == constrain_decreasing))
										//							throw new RuntimeException(Util.sprintf("Unexpected changepoint %s:%d:%c",chr,position,strand));

										// make sure that the binning did not result in a cp outside of the queried segment
										if(position >= start && position <= end){
											//								 if(segment_mle[i]>segment_mle[i+1]){
											//									 left.add(new AnnotatedRegion("exon", exon.chr, exon.start, position, exon.strand));
											//									 bw.write("iso", exon.chr, exon.start, position, exon.strand);
											//								 }
											//								 else{
											//									 right.add(new AnnotatedRegion("exon", exon.chr, exon.start, position, exon.strand));
											//									 bw.write("iso", exon.chr, position, exon.end, exon.strand);
											//								 }
											System.out.printf("%s %s %s\t%d\t%d\t%.2f\t%.2f\n", seg, segment_mle[i]>segment_mle[i+1]?"decrease":"increase",seg.chr+":"+position, Math.abs(seg.start-position), Math.abs(seg.end-position), segment_mle[i],segment_mle[i+1]);
											Map<String,Object> attributes = new HashMap<String, Object>();
											attributes.put("before_mle", segment_mle[i]);
											attributes.put("after_mle", segment_mle[i+1]);
											attributes.put("type","changepoint");

											changePoints.add(chr, position, position, strand, attributes);
											// write each changepoint
											AnnotatedRegion exon = IntervalTools.makeRegion(e.chr, e.get3Prime(), position, e.isNegativeStrand());
											gw.write(chr, exon.start, exon.end, strand);
											//							gw.write("exon",chr, position, position, strand, AnnotatedRegion.GTFAttributeString(attributes));
										}
									}
								}
							}
							
							
//							int[] y = BAMTools.binnedMaxEndCoverage(sfr, e.chr, seg.start, seg.size, 1, false, Strandedness.reverse_forward, e.isNegativeStrand());
//							double[] Y = new double[y.length];
//							for(int i=0; i<y.length; i++){
//								Y[i] = y[i];
//							}
//							IntegralTable integral = new IntegralTable(Y);
//							for(int i=0; i<y.length-1; i++){
//								double c1 = integral.sum(0, i);
//								double c2 = integral.sum(i+1, y.length-1);
//								
//								double P = (i+1.)/y.length;
//								double pdif = 1.0 - Beta.regularizedBeta(P, c2 + 1.0, c1);
//								if(pdif < alpha){
//									System.out.printf("%s\t%e\t%s:%d\n", seg,pdif, seg.chr, seg.start+i);	
//								}
//							}
						}
					}
				}
			}
		}
		gw.close();
	}

	public static void identifyComposite() throws FileNotFoundException{
		TranscriptIterator tiu = new TranscriptIterator("/home/sol/lailab/sol/GSE41637/300n/SRR594393.unspliced.gtf");
		TranscriptIterator ti = new TranscriptIterator("/home/sol/lailab/sol/GSE41637/300n/SRR594393.coverage.gtf");

		StrandedGenomicIntervalSet unspliced = IntervalTools.buildStrandedIntervalSet(tiu);
		StrandedGenomicIntervalTree<Map<String,Object>> isoscm = IntervalTools.buildRegionsTree(ti, false, true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> internal = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> spliced_exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> splice_jnct = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> splice_segs = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : isoscm){
			if("exon".equals(r.getAttribute("annotation"))){
				spliced_exons.add(r.chr, r.start, r.end, r.strand, r.attributes);
				if("internal_exon".equals(r.getAttribute("type"))){
					internal.add(r.chr, r.start, r.end, r.strand, r.attributes);
				}
			}
			if("splice_jnct".equals(r.getAttribute("annotation"))){
				splice_jnct.add(r.chr, r.start, r.end, r.strand, r.attributes);
				SegmentFragmenter.add(splice_segs, r.chr, r.start, r.end, r.strand);
			}
		}

		StrandedGenomicIntervalTree<Map<String,Object>> internal3p = IntervalTools.buildTerminiTree(internal, false, true, false);

		SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/GSE41637/mapped/indexed/SRR594393.bam"));

		//		for(AnnotatedRegion e : internal3p.overlappingRegions("15", 83115385, 83119232, '+')){
		for(AnnotatedRegion e : internal3p){
			List<AnnotatedRegion> splices = IntervalTools.BoundaryIntervals(splice_jnct, e.chr, IntervalTools.offsetPosition(e.get3Prime(),1,e.isNegativeStrand(),false), e.strand, true);
			List<AnnotatedRegion> segs = IntervalTools.BoundaryIntervals(splice_segs, e.chr, IntervalTools.offsetPosition(e.get3Prime(),1,e.isNegativeStrand(),false), e.strand, true);
			if(splices.size()>0 && segs.size()>0){
				Collections.sort(splices,new Comparator<AnnotatedRegion>() {
					public int compare(AnnotatedRegion o1, AnnotatedRegion o2) {
						return o1.size-o2.size;
					}
				});



				AnnotatedRegion shortest = splices.get(0);
				AnnotatedRegion seg = segs.get(0);
				List<AnnotatedRegion> un = IntervalTools.OverlappingIntervals(unspliced, e.chr, seg.start, seg.end, seg.strand);
				if(IntervalTools.OverlappingIntervals(spliced_exons, e.chr, seg.start, seg.end, seg.strand).size()==0 && un.size()==1){
					AnnotatedRegion us = un.get(0); 
					//				System.out.printf("segmenting %s\t%s %s\n", shortest,segs, us);

					String chr = seg.chr;
					int start = seg.start;
					int end = seg.end;
					int binSize=20;
					char strand = seg.strand;
					boolean isNegativeStrand=seg.isNegativeStrand();
					int maxBins=2000;
					Strandedness strandedness = Strandedness.reverse_forward;
					double alpha_0=1; double beta_0=1;
					int minCP=1;
					double nb_r=5;
					int r=10;
					double p=.99;
					//			boolean constrain_decreasing=false;
					double min_fold=.333;

					int l = end-start+1;

					int nBins = (l+binSize-1)/binSize;

					//			char strand = isNegativeStrand ? '-' : '+';


					if(nBins>maxBins){
						nBins=maxBins;
						binSize = (l+maxBins-1)/maxBins;
						System.out.printf("for segment %s:%d-%d (%c) adjusting binSize to %d\n", chr,start,end,strand,binSize);
					}

					// identify change points
					// cluster them, selecting the most likely from each cluster
					// write change points to GTF

					StrandedGenomicIntervalTree<Map<String,Object>> changePoints = new StrandedGenomicIntervalTree<Map<String,Object>>();



					// calculate the coverage in this chunk
					int[] cov = BAMTools.binnedMaxEndCoverage(sfr, chr, start, nBins, binSize, false, true, strandedness, isNegativeStrand);
					double[] y = new double[nBins];

					for(int i=0; i<nBins; i++){
						y[i] = cov[i];
					}

					//		SegmentationResult segmentation = TestConstraint.constrained_viterbi_segmentation(y, alpha_0, beta_0, nb_r, r, p, constrain_decreasing);
					SegmentationResult segmentation = IdentifyConstrainedChangePoints.doubly_constrained_segmentation(y, alpha_0, beta_0, nb_r, r, p, !seg.isNegativeStrand(), min_fold);

					if(segmentation.n_segments>1){
						//				int[][] segments = segmentation.segments();
						double[] segment_mle = segmentation.segment_mle;
						int[] changepoints = segmentation.change_point;

						for(int i=0 ; i<changepoints.length; i++){

							int changepoint_i = changepoints[i];
							// if it's an increase in coverage, the annotation should be shifted over by one bin
							if(segment_mle[i]<segment_mle[i+1]){
								changepoint_i++;					
							}

							if(changepoints[i] > minCP && changepoints[i] < nBins - minCP - 1){

								int position = start + (changepoint_i*binSize) + (binSize/2);

								//						if(!(segment_mle[i]>segment_mle[i+1] == constrain_decreasing))
								//							throw new RuntimeException(Util.sprintf("Unexpected changepoint %s:%d:%c",chr,position,strand));

								// make sure that the binning did not result in a cp outside of the queried segment
								if(position >= start && position <= end){
									//								 if(segment_mle[i]>segment_mle[i+1]){
									//									 left.add(new AnnotatedRegion("exon", exon.chr, exon.start, position, exon.strand));
									//									 bw.write("iso", exon.chr, exon.start, position, exon.strand);
									//								 }
									//								 else{
									//									 right.add(new AnnotatedRegion("exon", exon.chr, exon.start, position, exon.strand));
									//									 bw.write("iso", exon.chr, position, exon.end, exon.strand);
									//								 }
									System.out.printf("%s %s %s\t%d\t%.2f\t%.2f\n", seg, segment_mle[i]>segment_mle[i+1]?"decrease":"increase",seg.chr+":"+position, IntervalTools.makeRegion(seg.chr, seg.get5Prime(), position, seg.isNegativeStrand()).size, segment_mle[i],segment_mle[i+1]);
									Map<String,Object> attributes = new HashMap<String, Object>();
									attributes.put("before_mle", segment_mle[i]);
									attributes.put("after_mle", segment_mle[i+1]);
									attributes.put("type","changepoint");

									changePoints.add(chr, position, position, strand, attributes);
									// write each changepoint 
									//							gw.write("exon",chr, position, position, strand, AnnotatedRegion.GTFAttributeString(attributes));
								}
							}
						}
					}
					//				System.out.printf("%s\t%s\t%s\n", exon,StringUtils.join(left,","),StringUtils.join(right,","));

				}
				//				int[] i = IntervalTools.offsetInterval(shortest.start, shortest.end, -1, -1, e.isNegativeStrand());
				//				if(IntervalTools.BoundaryIntervals(splice_jnct, e.chr, i[0], i[1], e.strand).size()>0)
				//					System.out.printf("%s\t%s\n", splices.get(0),splices.get(splices.size()-1));	
			}
		}

	}

	public static void filterChangePoints(SAMFileReader sfr, File spliced_exon_gtf, File intronic_exon_gtf, File acc_jnct_gtf, File changepoint_gtf, GTFWriter filtered_changepoint_writer, int sj_radius) throws FileNotFoundException{
		StrandedGenomicIntervalTree<Map<String, Object>> continuousSegments = new StrandedGenomicIntervalTree<Map<String,Object>>();

		TranscriptIterator bi = new TranscriptIterator(acc_jnct_gtf);
		StrandedGenomicIntervalTree<Map<String,Object>> sj = IntervalTools.buildRegionsTree(bi, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> sj5p = IntervalTools.buildTerminiTree(sj, true, true, false);
		StrandedGenomicIntervalTree<Map<String,Object>> sj3p = IntervalTools.buildTerminiTree(sj, false, true, false);

		for(File gtf : Util.list(spliced_exon_gtf, intronic_exon_gtf)){
			TranscriptIterator ti = new TranscriptIterator(gtf);
			for(AnnotatedRegion segment : ti){
				SegmentFragmenter.add(continuousSegments, segment.chr, segment.start, segment.end, segment.strand);		
			}
		}

		StrandedGenomicIntervalTree<Map<String, Object>> ends5p = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> ends3p = new StrandedGenomicIntervalTree<Map<String,Object>>();

		TranscriptIterator ti = new TranscriptIterator(changepoint_gtf);
		for(AnnotatedRegion cp : ti){

			double before_mle = Double.parseDouble((String) cp.getAttribute("before_mle"));
			double after_mle = Double.parseDouble((String) cp.getAttribute("after_mle"));

			boolean isEnd5p = before_mle<after_mle ^ cp.isNegativeStrand();

			if(isEnd5p && !ends5p.contains(cp.chr, cp.start, cp.end, cp.strand) && !sj3p.overlappingRegions(cp.chr, cp.start-sj_radius, cp.end+sj_radius, cp.strand).iterator().hasNext()){
				ends5p.add(cp.chr, cp.start, cp.end, cp.strand, cp.attributes);
			}
			if(!isEnd5p && !ends3p.contains(cp.chr, cp.start, cp.end, cp.strand) && !sj5p.overlappingRegions(cp.chr, cp.start-sj_radius, cp.end+sj_radius, cp.strand).iterator().hasNext()){
				ends3p.add(cp.chr, cp.start, cp.end, cp.strand, cp.attributes);
			}
		}

		Comparator<AnnotatedRegion> descending_order = new Comparator<AnnotatedRegion>() {
			public int compare(AnnotatedRegion arg0, AnnotatedRegion arg1) {
				return arg0.start < arg1.start ? 1 : arg0.start > arg1.start ? -1 : 0;
			}
		};

		Comparator<AnnotatedRegion> ascending_order = new Comparator<AnnotatedRegion>() {
			public int compare(AnnotatedRegion arg0, AnnotatedRegion arg1) {
				return arg0.start < arg1.start ? -1 : arg0.start > arg1.start ? 1 : 0;
			}
		};

		for(AnnotatedRegion segment : continuousSegments){
			//			System.out.println(segment);

			{
				Comparator<AnnotatedRegion> c = segment.isNegativeStrand() ? ascending_order : descending_order;

				String attribute = segment.isNegativeStrand() ? "after_mle" : "before_mle";

				List<AnnotatedRegion> to_delete = new ArrayList<AnnotatedRegion>();
				List<AnnotatedRegion> overlapping3p = IntervalTools.OverlappingIntervals(ends3p, segment.chr, segment.start, segment.end, segment.strand);
				Collections.sort(overlapping3p, c);
				double max_mle = 0;
				for(AnnotatedRegion end : overlapping3p){
					double mle = Double.parseDouble((String) end.getAttribute(attribute));
					if(mle > max_mle){
						max_mle = Math.max(mle, max_mle);
						//					System.out.printf("end %s %s\n", end, end.attributes);
					}
					else{
						to_delete.add(end);
					}
				}

				for(AnnotatedRegion end : to_delete){
					ends3p.remove(end.chr, end.start, end.end, end.strand);
				}
			}

			{
				Comparator<AnnotatedRegion> c = segment.isNegativeStrand() ? descending_order : ascending_order;

				String attribute = segment.isNegativeStrand() ? "before_mle" : "after_mle";

				List<AnnotatedRegion> to_delete = new ArrayList<AnnotatedRegion>();
				List<AnnotatedRegion> overlapping5p = IntervalTools.OverlappingIntervals(ends5p, segment.chr, segment.start, segment.end, segment.strand);
				Collections.sort(overlapping5p, c);
				double max_mle = 0;
				for(AnnotatedRegion end : overlapping5p){
					double mle = Double.parseDouble((String) end.getAttribute(attribute));
					if(mle > max_mle){
						max_mle = Math.max(mle, max_mle);
					}
					else{
						to_delete.add(end);
					}
				}

				for(AnnotatedRegion end : to_delete){
					ends5p.remove(end.chr, end.start, end.end, end.strand);
				}
			}

			//				System.out.println(overlapping3p);
			//				System.out.println(to_delete);
		}

		for(AnnotatedRegion end : ends3p){
			filtered_changepoint_writer.write("exon", end.chr, end.start, end.start, end.strand, AnnotatedRegion.GTFAttributeString(end.attributes));
		}
		for(AnnotatedRegion end : ends5p){
			filtered_changepoint_writer.write("exon", end.chr, end.start, end.start, end.strand, AnnotatedRegion.GTFAttributeString(end.attributes));
		}

	}

	public static void enumerateChangepointExons(File spliced_exon_gtf, File intronic_exon_gtf, File changepoint_gtf, GTFWriter changepoint_exon_writer) throws FileNotFoundException {
		StrandedGenomicIntervalTree<Map<String, Object>> continuousSegments = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> written = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(File gtf : Util.list(spliced_exon_gtf, intronic_exon_gtf)){
			TranscriptIterator ti = new TranscriptIterator(gtf);
			for(AnnotatedRegion segment : ti){
				continuousSegments.add(segment.chr, segment.start, segment.end, segment.strand, segment.attributes);	
			}
		}

		TranscriptIterator ti = new TranscriptIterator(changepoint_gtf);
		for(AnnotatedRegion cp : ti){

			double before_mle = Double.parseDouble((String) cp.getAttribute("before_mle"));
			double after_mle = Double.parseDouble((String) cp.getAttribute("after_mle"));

			boolean isEnd5p = before_mle<after_mle ^ cp.isNegativeStrand();

			for(AnnotatedRegion segment : continuousSegments.overlappingRegions(cp.chr, cp.start, cp.end, cp.strand)){
				String type = (String) segment.getAttribute("type");
				if(type.equals("5p_exon") && isEnd5p){
					// write the 5p exon
					if(!written.contains(segment.chr, Math.min(segment.get3Prime(),cp.start), Math.max(segment.get3Prime(),cp.start), cp.strand)){
						written.add(segment.chr, Math.min(segment.get3Prime(),cp.start), Math.max(segment.get3Prime(),cp.start), cp.strand);
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("type", "5p_exon");
						changepoint_exon_writer.write("5p_exon", segment.chr, Math.min(segment.get3Prime(),cp.start), Math.max(segment.get3Prime(),cp.start), cp.strand, AnnotatedRegion.GTFAttributeString(attributes));
					}
				}
				if(type.equals("3p_exon") && !isEnd5p){
					// write the 3p exon
					if(!written.contains(segment.chr, Math.min(segment.get5Prime(),cp.start), Math.max(segment.get5Prime(),cp.start), cp.strand)){
						written.add(segment.chr, Math.min(segment.get5Prime(),cp.start), Math.max(segment.get5Prime(),cp.start), cp.strand);
						Map<String,Object> attributes = new HashMap<String, Object>();
						attributes.put("type", "3p_exon");
						changepoint_exon_writer.write("3p_exon", segment.chr, Math.min(segment.get5Prime(),cp.start), Math.max(segment.get5Prime(),cp.start), cp.strand, AnnotatedRegion.GTFAttributeString(attributes));
					}
				}
				if(type.equals("internal_exon")){
					// if it's the 5p or the 3p, whichever is appropriate
					if(isEnd5p){
						if(!written.contains(segment.chr, Math.min(segment.get3Prime(),cp.start), Math.max(segment.get3Prime(),cp.start), cp.strand)){
							written.add(segment.chr, Math.min(segment.get3Prime(),cp.start), Math.max(segment.get3Prime(),cp.start), cp.strand);
							Map<String,Object> attributes = new HashMap<String, Object>();
							attributes.put("type", "5p_exon");
							changepoint_exon_writer.write("5p_exon", segment.chr, Math.min(segment.get3Prime(),cp.start), Math.max(segment.get3Prime(),cp.start), cp.strand, AnnotatedRegion.GTFAttributeString(attributes));
						}
					}
					else{
						if(!written.contains(segment.chr, Math.min(segment.get5Prime(),cp.start), Math.max(segment.get5Prime(),cp.start), cp.strand)){
							written.add(segment.chr, Math.min(segment.get5Prime(),cp.start), Math.max(segment.get5Prime(),cp.start), cp.strand);
							Map<String,Object> attributes = new HashMap<String, Object>();
							attributes.put("type", "3p_exon");
							changepoint_exon_writer.write("3p_exon", segment.chr, Math.min(segment.get5Prime(),cp.start), Math.max(segment.get5Prime(),cp.start), cp.strand, AnnotatedRegion.GTFAttributeString(attributes));
						}
					}
				}
			}	
		}
	}

	public static void main(String[] args) throws IOException {
//		if(true){
//			File dir = new File("/home/sol/data/GSE41637/isoscm/SRR594393");
//			String base = "SRR594393";
//			BBFileReader pos_cov = new BBFileReader("/home/sol/lailab/sol/GSE41637/bigwig/bigwig/SRR594393.plus.bw");
//			BBFileReader neg_cov = new BBFileReader("/home/sol/lailab/sol/GSE41637/bigwig/bigwig/SRR594393.minus.bw");
//			
//			File assembly_gtf = new File("/home/sol/data/GSE41637/isoscm/SRR594393/SRR594393.gtf");
//			Strandedness strandedness = Strandedness.reverse_forward;
//			System.out.println("calculating coverage");
//			File coverage_gtf = new File(Util.sprintf("%s/%s.coverage.gtf",dir, base));
//		
//			GTFWriter coverage_gw = new GTFWriter(IO.bufferedPrintstream(coverage_gtf));
//			IdentifyChangePoints.calculate_spliced_coverage(pos_cov,neg_cov, assembly_gtf, strandedness, coverage_gw);
//			coverage_gw.close();
//			System.exit(0);
//		}
		if(false){
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			GTFWriter gw = new GTFWriter("/dev/stdout");
			File unspliced_gtf = new File("/home/sol/lailab/sol/GSE41637/300n/SRR594393.unspliced.gtf");
			File spliced_gtf = new File("/home/sol/lailab/sol/GSE41637/300n/tmp/SRR594393.exon.gtf");
			int minLength=100;
			int maxBins=2000;
			int binSize=5;
			int minCP=0;
			Strandedness strandedness = Strandedness.reverse_forward;
			double alpha_0=1; double beta_0=1;
			int nb_r=1;
			int r=10;
			double p=.99;
			double min_fold=.33;
			int min_terminal=-100;
			identifyNegativeBinomialChangePointsInUnsplicedSegments(sfr, unspliced_gtf, spliced_gtf, gw, minLength, maxBins, binSize, minCP, strandedness, alpha_0, beta_0, nb_r, r, p, min_fold, min_terminal);
		}
		if(true){
			//			identifyComposite();
			identifyComposite2();
		}
		
		if(false){
			SAMFileReader sfr = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			File unspliced_gtf = new File("/home/sol/lailab/sol/GSE41637/nb/gtf/SRR594393.remainder.gtf");
			File assembly_gtf = new File("/home/sol/lailab/sol/GSE41637/nb/gtf/SRR594393.gtf");
			GTFWriter fcw = new GTFWriter("/dev/stdout");

			//			calculate_spliced_coverage(sfr, assembly_gtf, Strandedness.reverse_forward, fcw);
			calculate_coverage(sfr, unspliced_gtf, Strandedness.reverse_forward, fcw);
		}
		if(false){
			SAMFileReader sfr = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			File spliced_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.1.exon.gtf");
			File intronic_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.1.intronic.gtf");
			File changepoint_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.1.cp_exon.gtf");
			GTFWriter fcw = new GTFWriter("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.1.cp_exon.selected.gtf");

			double m_0=0;
			double V_0=1; 
			double a_0=1; 
			double b_0=1;

			//			int r = 50; 
			int r = 25*10; 
			double p = .95;

//			selectChangePoints(sfr, spliced_exon_gtf, intronic_exon_gtf, changepoint_exon_gtf, 3, Strandedness.reverse_forward, m_0, V_0, a_0, b_0, r, p, fcw);
			fcw.close();
		}

		if(false){
			SAMFileReader sfr1 = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594393.bam"));
			File spliced_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.exon.gtf");
			File sj_bed = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.sj.bed");
			File intronic_exon_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.intronic.gtf");
			File changepoint_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.cp.gtf");
			File filtered_changepoint_gtf = new File("/home/sol/workspace/IsoSCM/tests/gtf/tmp/SRR594393.25.cp.filtered.gtf");
			GTFWriter fcw = new GTFWriter(new PrintStream(filtered_changepoint_gtf));
			int minLength = 50;
			filterChangePoints(sfr1, spliced_exon_gtf, intronic_exon_gtf, sj_bed, changepoint_gtf, fcw, minLength);
			fcw.close();
		}

		if(false){
			SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/GSE41637/isoscm/bam/enah.brain.bam"));
			double alpha = 1.5;
			double beta = .05;
			int r = 50; 
			//		double p = .95;
			double p = .99;

			int nBins=100;
			int binSize=25;
			int binShift=8;
			//		int binShift=123;
			int binClusterWidth = 5;
			int minCP=35;
			GTFWriter gw = new GTFWriter("/home/sol/lailab/sol/GSE41637/isoscm/gtf/tmp/enah.brain.cp.gtf");
//			identifyPoissonChangePoints(sfr, "16", 29229680, 29239024, nBins, binSize, binShift, minCP, binClusterWidth, Strandedness.reverse_forward, true, alpha, beta, r, p, gw);
			gw.close();
		}
	}
}
