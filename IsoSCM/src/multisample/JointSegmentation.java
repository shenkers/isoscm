package multisample;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import scm.IdentifyChangePoints;
import tools.AnnotatedRegion;
import tools.GTFTools.GTFWriter;
import tools.IntervalTools;
import tools.ParseGTF.TranscriptIterator;
import tools.StrandedGenomicIntervalSet;
import tools.StrandedGenomicIntervalTree;
import tools.Strandedness;
import util.IO;
import util.Util;
import util.Util.Criteria;
import util.Util.Function;
import changepoint.ChangePoint;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class JointSegmentation {
	private static final Logger logger = LogManager.getLogger();

	public static void performJointSegmentationUsingReference(String ensemblf, String id1, String id2, String exons1, String exons2, String bam1, String bam2, String outFile) throws FileNotFoundException {

		StrandedGenomicIntervalTree<Map<String,Object>> ensembl = IntervalTools.buildRegionsTree(new TranscriptIterator(new File(ensemblf)), true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> isoscm1 = IntervalTools.buildRegionsTree(new TranscriptIterator(new File(exons1)), true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> isoscm2 = IntervalTools.buildRegionsTree(new TranscriptIterator(new File(exons2)), true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> t5p = IntervalTools.buildTerminiTree(ensembl, true, true, true);

		int maxBins = 20000;
		int binSize=10;
		int minCP=0;
		Strandedness strandedness = Strandedness.reverse_forward;
		double alpha_0=1;
		double beta_0=1; 
		int nb_r=1; 
		int r=10;
		double p=.99;

		double min_fold=.6;

		String[] ids = new String[]{
				id1,
				id2
		};
		SamReader[] sfrs = new SamReader[]{
				SamReaderFactory.makeDefault().open(new File(bam1)),
				SamReaderFactory.makeDefault().open(new File(bam2)),
		}; 

		GTFWriter gw = new GTFWriter(outFile);
		int id=0;
		for(AnnotatedRegion region : t5p){
			//				if(!region.chr.equals("10") || !IntervalTools.isContained(region.start, region.end, 95541467, 95545255))
			//					continue;

			List<AnnotatedRegion> iso_ends1 = IntervalTools.BoundaryIntervals(isoscm1, region.chr, region.get5Prime(), region.strand, true);
			List<AnnotatedRegion> iso_ends2 = IntervalTools.BoundaryIntervals(isoscm2, region.chr, region.get5Prime(), region.strand, true);

			if(iso_ends1.size()>0 && iso_ends2.size()>0){
				//					System.out.println(iso_ends1);
				//					System.out.println(iso_ends2);
				StrandedGenomicIntervalSet union = IntervalTools.buildStrandedIntervalSet(iso_ends1);
				IntervalTools.addRegions(union, iso_ends2);
				for(AnnotatedRegion u : union){
					Map<String,Object> attributes = new HashMap<String, Object>();
					String ID = Util.sprintf("%09d", id);
					attributes.put("id", ID);
					id++;
					gw.write("locus", u.chr, u.start, u.end, u.strand, AnnotatedRegion.GTFAttributeString(attributes));
					//System.out.println(u);
					//				else if("3p_exon".equals(type)){
					boolean constrained_decreasing = !u.isNegativeStrand();
					StrandedGenomicIntervalTree<Map<String, Object>> msr = null;//IdentifyChangePoints.identifyConstrainedNegativeBinomialPoints(ids, sfrs, u.chr, u.start, u.end, maxBins, binSize, minCP, strandedness, u.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, constrained_decreasing, min_fold);

					// write the union exon
					for(AnnotatedRegion cp : msr){
						cp.addAttribute("s1", id1);
						cp.addAttribute("s2", id2);
						cp.addAttribute("before_mle", StringUtils.join(Util.list(cp.getAttribute("before_mle")),","));
						cp.addAttribute("after_mle", StringUtils.join(Util.list(cp.getAttribute("after_mle")),","));
						cp.addAttribute("id", ID);
						gw.write("changepoint", u.chr, cp.start, cp.end, u.strand, AnnotatedRegion.GTFAttributeString(cp.attributes));
						// write all change points, we'll filter out later when we decide on parameters...
					}

				}
			}
		}

		gw.close();

	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> get3pTerminalExons(StrandedGenomicIntervalTree<Map<String,Object>> spliced_exons){
		StrandedGenomicIntervalTree<Map<String,Object>> terminal_exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		for(AnnotatedRegion exon : spliced_exons){
			if("3p_exon".equals(exon.getAttribute("type"))){
				terminal_exons.add(exon);
			}
		}
		return terminal_exons;
	}

	public static StrandedGenomicIntervalTree<Map<String, Object>> getMaxCommonExons(StrandedGenomicIntervalTree<Map<String,Object>> exons1, StrandedGenomicIntervalTree<Map<String,Object>> exons2, Strandedness s1, Strandedness s2){
		// get the terminal 3p exons
		StrandedGenomicIntervalTree<Map<String, Object>> t3p1 = get3pTerminalExons(exons1);
		StrandedGenomicIntervalTree<Map<String, Object>> t3p2 = get3pTerminalExons(exons2);
		// get their 5p ends
		StrandedGenomicIntervalTree<Map<String, Object>> t5p = IntervalTools.buildTerminiTree(t3p1, true, true, false);
		IntervalTools.addRegions(t5p, IntervalTools.buildTerminiTree(t3p2, true, true, false), true, false);

		StrandedGenomicIntervalTree<Map<String,Object>> maxCommonExon = new StrandedGenomicIntervalTree<Map<String,Object>>();
		for(AnnotatedRegion e5p : t5p){
			// check that both samples have this terminal exon
			List<AnnotatedRegion> e3ps1 = IntervalTools.BoundaryIntervals(t3p1, e5p.chr, e5p.start, e5p.strand, true);
			List<AnnotatedRegion> e3ps2 = IntervalTools.BoundaryIntervals(t3p2, e5p.chr, e5p.start, e5p.strand, true);
			if(e3ps1.size()>0 && e3ps2.size()>0){
				AnnotatedRegion e3p1 = e3ps1.get(0);
				AnnotatedRegion e3p2 = e3ps2.get(0);
				AnnotatedRegion maxExon = new AnnotatedRegion("", e5p.chr, Math.min(e3p1.start,e3p2.start), Math.max(e3p1.end,e3p2.end), e5p.strand);
				//make sure the maxExon overlaps at most 1 exon in each sample, so that we aren't merging close genes
				if(IntervalTools.OverlappingIntervals(exons1, maxExon.chr, maxExon.start, maxExon.end, maxExon.strand).size()==1&& IntervalTools.OverlappingIntervals(exons2, maxExon.chr, maxExon.start, maxExon.end, maxExon.strand).size()==1){
					// if either of the datasets is unstranded
					if(s1==Strandedness.unstranded||s2==Strandedness.unstranded){
						// make sure the max exon doesn't overlap anything on the opposite strand
						if(!(exons1.overlappingRegions(maxExon.chr, maxExon.start, maxExon.end, IntervalTools.opposite(maxExon.strand)).iterator().hasNext()||exons2.overlappingRegions(maxExon.chr, maxExon.start, maxExon.end, IntervalTools.opposite(maxExon.strand)).iterator().hasNext())){
							maxCommonExon.add(maxExon);						
						}
					}
					else{
						maxCommonExon.add(maxExon);
					}
				}
			}
		}

		return maxCommonExon;
	}

	public static void performJointSegmentation(String id1, String id2, File spliced_exons1, File spliced_exons2, File bam1, File bam2, File table, File gtf, Strandedness strandedness1, Strandedness strandedness2, int maxBins, int binSize, int minCP, double alpha_0, double beta_0, int nb_r, int r, double p, double min_fold) throws FileNotFoundException {

		StrandedGenomicIntervalTree<Map<String,Object>> exons1 = IntervalTools.buildRegionsTree(new TranscriptIterator(spliced_exons1), true, true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> exons2 = IntervalTools.buildRegionsTree(new TranscriptIterator(spliced_exons2), true, true, true);

		logger.info("identifying common terminal exons...");
		StrandedGenomicIntervalTree<Map<String, Object>> maxCommonExons = getMaxCommonExons(exons1, exons2, strandedness1, strandedness2);
		StrandedGenomicIntervalTree<Map<String, Object>> t5p = IntervalTools.buildTerminiTree(maxCommonExons, true, true, false);

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		String[] ids = new String[]{
				id1,
				id2
		};
		SamReader[] sfrs = new SamReader[]{
				SamReaderFactory.makeDefault().open(bam1),
				SamReaderFactory.makeDefault().open(bam2),
		};
		Strandedness[] strandednesses = new Strandedness[]{
				strandedness1,
				strandedness2,
		};
		
		logger.info("Performing joint segmentation...");
		GTFWriter gw = new GTFWriter(IO.bufferedPrintstream(gtf));
		PrintStream tabular = IO.bufferedPrintstream(table);
		tabular.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "samples","locus_id","changepoint","confidence","log_odds","upstream_segment","downstream_segment","locus","strand","upstream_cov","downstream_cov","site_usage","differential_usage");
		int changepoint_id=0;
		for(AnnotatedRegion region : t5p){
			
			AnnotatedRegion u = IntervalTools.BoundaryIntervals(maxCommonExons, region.chr, region.get5Prime(), region.strand, true).iterator().next();

			Map<String,Object> attributes = new HashMap<String, Object>();
			String locus_ID = Util.sprintf("%09d", changepoint_id);
			attributes.put("id", locus_ID);
			changepoint_id++;
			gw.write("locus", u.chr, u.start, u.end, u.strand, AnnotatedRegion.GTFAttributeString(attributes));
		
			boolean constrained_decreasing = !u.isNegativeStrand();
			
			List<ChangePoint> changepoints = IdentifyChangePoints.identifyConstrainedNegativeBinomialPoints(ids,sfrs,strandednesses, u.chr, u.start, u.end, maxBins, binSize, minCP, u.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, constrained_decreasing, min_fold);

			// write the union exon
			for(ChangePoint cp : changepoints){
				AnnotatedRegion changepoint = new AnnotatedRegion("changepoint", cp.pos.chr, cp.pos.start, cp.pos.end, cp.pos.strand, new HashMap<String, Object>());
				changepoint.addAttribute("s1", id1);
				changepoint.addAttribute("s2", id2);
				changepoint.addAttribute("before_mle", StringUtils.join(Util.list(changepoint.getAttribute("before_mle")),","));
				changepoint.addAttribute("after_mle", StringUtils.join(Util.list(changepoint.getAttribute("after_mle")),","));
				changepoint.addAttribute("cov_upstream", StringUtils.join(Util.list(cp.cov_upstream),","));
				changepoint.addAttribute("cov_downstream", StringUtils.join(Util.list(cp.cov_downstream),","));
				changepoint.addAttribute("confidence", Util.sprintf("%.3e",cp.confidence));
				changepoint.addAttribute("log_odds", Util.sprintf("%.3e",cp.log_odds));
				changepoint.addAttribute("id", locus_ID);
				gw.write("changepoint", u.chr, changepoint.start, changepoint.end, u.strand, AnnotatedRegion.GTFAttributeString(changepoint.attributes));


				double[] usage = new double[sfrs.length];
				for (int i = 0; i < sfrs.length; i++) {
					usage[i] =  1-Math.min((cp.cov_downstream[i]/cp.cov_upstream[i]),1);
				}

				// write all change points, we'll filter out later when we decide on parameters...
				tabular.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", StringUtils.join(ids,","),locus_ID,changepoint,Util.sprintf("%.3e",cp.confidence), Util.sprintf("%.3e",cp.log_odds),cp.upstream_region,cp.downstream_region,u,cp.pos.strand,StringUtils.join(Util.list(cp.cov_upstream),","),StringUtils.join(Util.list(cp.cov_downstream),","),StringUtils.join(Util.list(usage),","),usage[0]-usage[1]);

			}

		}

		gw.close();
		tabular.close();
	}

	public static void identifyBypassedRegions(File f) throws IOException{

		TranscriptIterator ti = new TranscriptIterator(f);
		StrandedGenomicIntervalTree<Map<String,Object>> t = IntervalTools.buildRegionsTree(ti, false, true, true);

		System.out.printf("location\tpcov1\tpcov2\tbypass1-bypass2\tbypass1\tbypass2\n");
		for(AnnotatedRegion region : t){
			if(region.getAttribute("annotation").equals("locus")){

				for(AnnotatedRegion cp : IntervalTools.OverlappingIntervals(t, region.chr, region.start, region.end, region.strand)){

					if(cp.getAttribute("annotation").equals("locus"))
						continue;
					if(Math.abs(cp.start-region.get3Prime())<300){
						continue;
					}



					//						System.out.println(cp.toAttributeString());
					Function<String,Double> cast = new Function<String,Double>(){
						public Double evaluate(String s) {
							return Double.parseDouble(s);
						}
					};
					List<Double> before = Util.evaluate(Util.list(((String) cp.getAttribute("before_mle")).split(",")), cast);
					List<Double> after = Util.evaluate(Util.list(((String) cp.getAttribute("after_mle")).split(",")), cast);

					double prox0 = cp.isNegativeStrand() ? after.get(0): before.get(0);
					double prox1 = cp.isNegativeStrand() ? after.get(1): before.get(1);

					double r0 = cp.isNegativeStrand() ? before.get(0)/after.get(0) : after.get(0)/before.get(0);
					double r1 = cp.isNegativeStrand() ? before.get(1)/after.get(1) : after.get(1)/before.get(1);
					r0 = r0 > 1 ? 1 : r0;
					r1 = r1 > 1 ? 1 : r1;

					System.out.printf("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", cp, prox0, prox1, r0-r1, r0, r1);

				}
			}
		}

	}
}
