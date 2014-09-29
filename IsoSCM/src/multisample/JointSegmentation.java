package multisample;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.apache.commons.lang.StringUtils;

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
import util.Util.Function;
import changepoint.ChangePoint;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

public class JointSegmentation {

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
		SAMFileReader[] sfrs = new SAMFileReader[]{
				new SAMFileReader(new File(bam1)),
				new SAMFileReader(new File(bam2)),
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
	
	public static void performJointSegmentation(String id1, String id2, String exons1, String exons2, String bam1, String bam2, String outTabularFile, String outGtfFile, Strandedness strandedness) throws FileNotFoundException {

		StrandedGenomicIntervalTree<Map<String,Object>> isoscm1 = IntervalTools.buildRegionsTree(new TranscriptIterator(new File(exons1)), true, true, true);
		StrandedGenomicIntervalTree<Map<String,Object>> isoscm2 = IntervalTools.buildRegionsTree(new TranscriptIterator(new File(exons2)), true, true, true);

		StrandedGenomicIntervalTree<Map<String,Object>> t5p = new StrandedGenomicIntervalTree<Map<String,Object>>();
		for(StrandedGenomicIntervalTree<Map<String,Object>> isoscm : Util.list(isoscm1,isoscm2)){
		for(AnnotatedRegion r : isoscm){
//			System.out.println(r.toAttributeString());
			if("3p_exon".equals(r.getAttribute("type"))){
				if(!t5p.contains(r.chr, r.get5Prime(), r.get5Prime(), r.strand))
					t5p.add(r.chr, r.get5Prime(), r.get5Prime(), r.strand);
			}
		}
		}
		
		int maxBins = 20000;
		int binSize=10;
		int minCP=0;
		double alpha_0=1;
		double beta_0=1; 
		int nb_r=1; 
		int r=10;
		double p=.99;

		double min_fold=.6;

		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		String[] ids = new String[]{
				id1,
				id2
		};
		SAMFileReader[] sfrs = new SAMFileReader[]{
				new SAMFileReader(new File(bam1)),
				new SAMFileReader(new File(bam2)),
		}; 

		GTFWriter gw = new GTFWriter(outGtfFile);
		PrintStream tabular = IO.bufferedPrintstream(outTabularFile);
		tabular.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "samples","locus_id","changepoint","upstream_segment","downstream_segment","locus","strand","upstream_cov","downstream_cov","site_usage","differential_usage");
		int changepoint_id=0;
		for(AnnotatedRegion region : t5p){
			//				if(!region.chr.equals("10") || !IntervalTools.isContained(region.start, region.end, 95541467, 95545255))
			//					continue;

			List<AnnotatedRegion> iso_ends1 = IntervalTools.BoundaryIntervals(isoscm1, region.chr, region.get5Prime(), region.strand, true);
			List<AnnotatedRegion> iso_ends2 = IntervalTools.BoundaryIntervals(isoscm2, region.chr, region.get5Prime(), region.strand, true);

			for(List<AnnotatedRegion> i : Util.list(iso_ends1, iso_ends2)){
				Iterator<AnnotatedRegion> it = i.iterator();
				while(it.hasNext()){
					AnnotatedRegion n = it.next();
					if(!"3p_exon".equals(n.getAttribute("type")))
						it.remove();
				}
			}
			if(iso_ends1.size()>0 && iso_ends2.size()>0){
				//					System.out.println(iso_ends1);
				//					System.out.println(iso_ends2);
				StrandedGenomicIntervalSet union = IntervalTools.buildStrandedIntervalSet(iso_ends1);
				IntervalTools.addRegions(union, iso_ends2);
				for(AnnotatedRegion u : union){
					Map<String,Object> attributes = new HashMap<String, Object>();
					String locus_ID = Util.sprintf("%09d", changepoint_id);
					attributes.put("id", locus_ID);
					changepoint_id++;
					gw.write("locus", u.chr, u.start, u.end, u.strand, AnnotatedRegion.GTFAttributeString(attributes));
					//System.out.println(u);
					//				else if("3p_exon".equals(type)){
					boolean constrained_decreasing = !u.isNegativeStrand();
//					2R:454973-454973
					 
					List<ChangePoint> changepoints = IdentifyChangePoints.identifyConstrainedNegativeBinomialPoints(ids,sfrs, u.chr, u.start, u.end, maxBins, binSize, minCP, strandedness, u.isNegativeStrand(), alpha_0, beta_0, nb_r, r, p, constrained_decreasing, min_fold);

					// write the union exon
					for(ChangePoint cp : changepoints){
						AnnotatedRegion changepoint = new AnnotatedRegion("changepoint", cp.pos.chr, cp.pos.start, cp.pos.end, cp.pos.strand, new HashMap<String, Object>());
						changepoint.addAttribute("s1", id1);
						changepoint.addAttribute("s2", id2);
						changepoint.addAttribute("before_mle", StringUtils.join(Util.list(changepoint.getAttribute("before_mle")),","));
						changepoint.addAttribute("after_mle", StringUtils.join(Util.list(changepoint.getAttribute("after_mle")),","));
						changepoint.addAttribute("cov_upstream", StringUtils.join(Util.list(cp.cov_upstream),","));
						changepoint.addAttribute("cov_downstream", StringUtils.join(Util.list(cp.cov_downstream),","));
						changepoint.addAttribute("id", locus_ID);
						gw.write("changepoint", u.chr, changepoint.start, changepoint.end, u.strand, AnnotatedRegion.GTFAttributeString(changepoint.attributes));
						
						
						double[] usage = new double[sfrs.length];
						for (int i = 0; i < sfrs.length; i++) {
							usage[i] = 1-(cp.cov_downstream[i]/cp.cov_upstream[i]);
						}
						
						// write all change points, we'll filter out later when we decide on parameters...
						tabular.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", StringUtils.join(ids,","),locus_ID,changepoint,cp.upstream_region,cp.downstream_region,u,cp.pos.strand,StringUtils.join(Util.list(cp.cov_upstream),","),StringUtils.join(Util.list(cp.cov_downstream),","),StringUtils.join(Util.list(usage),","),usage[0]-usage[1]);
						
					}

				}
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

	public static void main(String[] args) throws IOException {
		
		class MultiSegmentCommand{	
		
			@Parameter(names="-exons1", description="exons from sample 1")
			String exons1;

			@Parameter(names="-exons2", description="exons from sample 2")
			String exons2;

			@Parameter(names="-bam1", description="bam from sample 1")
			String bam1;

			@Parameter(names="-bam2", description="bam from sample 2")
			String bam2;

			@Parameter(names="-id1", description="id from sample 1")
			String id1;

			@Parameter(names="-id2", description="id from sample 2")
			String id2;

			@Parameter(names="-out", description="gtf to which output will be written")
			String out;

			@Parameter(names="-s", description="strandedness")
			String strandedness;
		}
		
		class ListCommand{	
			
			@Parameter(names="-pairfile", description="gtf from the pairwise analysis")
			String pairfile;

		}

		class HelpCommand{	

		}

		JCommander jc = new JCommander();
		MultiSegmentCommand compare = new MultiSegmentCommand();
		ListCommand list = new ListCommand();
		HelpCommand help = new HelpCommand();
		jc.addCommand("compare", compare);
		jc.addCommand("list", list);
		jc.addCommand("-h", help);
		jc.parse(args);

		if(jc.getParsedCommand()==null || jc.getParsedCommand().equals("-h")){
			jc.usage();
		}
		else if(jc.getParsedCommand().equals("compare")){
			performJointSegmentation(compare.id1, compare.id2, compare.exons1, compare.exons2, compare.bam1, compare.bam2, compare.out, null, Strandedness.valueOf(compare.strandedness));	
		}
		else if(jc.getParsedCommand().equals("list")){
			identifyBypassedRegions(new File(list.pairfile));	
		}


	}
}
