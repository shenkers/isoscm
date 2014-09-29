package executable;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import multisample.JointSegmentation;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.apache.commons.io.FileUtils;

import processing.ClusterExpressedSegments;
import processing.FindSpliceJunctions;
import processing.SlidingWindow;
import scm.IdentifyChangePoints;
import splicegraph.ExonSpliceGraph;
import tools.BEDTools.BEDWriter;
import tools.GTFTools.GTFWriter;
import tools.Strandedness;
import util.IO;
import util.Util;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.BooleanConverter;
import com.beust.jcommander.converters.DoubleConverter;
import com.beust.jcommander.converters.FileConverter;
import com.beust.jcommander.converters.IntegerConverter;

public class IsoSCM {

	public static void main(String[] args) throws IOException {
		class AssembleCommand{	
			@Parameter(names="-bam", description="bam file to be segmented", converter=FileConverter.class)
			File bam;

			@Parameter(names="-stringency", description="stringency for BAM format checking")
			String stringency;

			@Parameter(names="-t", description="the threshold above which a segment is considered expressed", converter=IntegerConverter.class)
			int threshold=1;

			@Parameter(names="-s", description="the strandedness, can be \"reverse_forward\" or \"unstranded\"")
			String strandedness;

			@Parameter(names="-dir", description="The output directory")
			String dir = "isoscm";

			@Parameter(names="-base", description="The output basename. This will be used as the prefix for output files.")
			String base;

			@Parameter(names="-insert_size_quantile", description="If the data is paired, will attempt to scaffold gaps between segments spanned by mates separated by this quantile")
			Double insert_size_quantile;

			@Parameter(names="-w", description="the width of window to be used by the SCM", converter=IntegerConverter.class)
			int w = 20;

			@Parameter(names="-segment_r", description="controls how long the fragments are", converter=IntegerConverter.class)
			int segment_r = 10;

			// can hide parameters by using hidden=true
			@Parameter(names="-segment_p", description="controls how long the fragments are", converter=DoubleConverter.class)
			double segment_p = 0.95;

			@Parameter(names="-nb_r", description="controls expected noisy-ness of the data", converter=IntegerConverter.class)
			int nb_r = 10;

			@Parameter(names="-merge_radius", description="gaps smaller than this distance will be merged", converter=IntegerConverter.class)
			int merge_radius = 100;

			//			@Parameter(names="-max_cps", description="maximum number of changepoints allowed in a continuous segment", converter=IntegerConverter.class)
			//			int max_cps;

			@Parameter(names="-internal", description="whether or not to identify internal change points", converter=BooleanConverter.class)
			boolean internal;

			@Parameter(names="-coverage", description="whether or not to calculate coverage of resulting models", converter=BooleanConverter.class, arity=1)
			boolean coverage=true;

			@Parameter(names="-min_fold", description="the minimum fold change between neighboring segments expressed as a ratio of low/high, acceptable values in the range [0.0-1.0]", converter=DoubleConverter.class)
			double min_fold = 0.5;

			@Parameter(names="-min_terminal", description="terminal segments are \"virtually\" extended by this amount before segmentation", converter=IntegerConverter.class)
			int min_terminal = 300;

			@Parameter(names="-jnct_alpha", description="The significance level for binomial test to accept splice junction", converter=DoubleConverter.class)
			double jnct_alpha = 0.05;

			@Parameter(names="-merge_segments", description="Optional:bed file of regions to merge across")
			String filled_gap_segments;

		}

		class SegmentCommand{	
			@Parameter(names="-bam", description="bam file to be segmented", converter=FileConverter.class)
			File bam;

			@Parameter(names="-s", description="the strandedness, can be \"reverse_forward\" or \"unstranded\"")
			String strandedness;

			@Parameter(names="-dir", description="The output directory")
			String dir = "isoscm";

			@Parameter(names="-base", description="The output basename. This will be used as the prefix for output files.")
			String base;

			@Parameter(names="-w", description="the width of window to be used by the SCM", converter=IntegerConverter.class)
			int w = 20;

			@Parameter(names="-segment_r", description="controls how long the fragments are", converter=IntegerConverter.class)
			int segment_r = 10;

			// can hide parameters by using hidden=true
			@Parameter(names="-segment_p", description="controls how long the fragments are", converter=DoubleConverter.class)
			double segment_p = 0.95;

			@Parameter(names="-nb_r", description="controls expected noisy-ness of the data", converter=IntegerConverter.class)
			int nb_r = 10;

			@Parameter(names="-merge_radius", description="gaps smaller than this distance will be merged", converter=IntegerConverter.class)
			int merge_radius = 100;

			//			@Parameter(names="-max_cps", description="maximum number of changepoints allowed in a continuous segment", converter=IntegerConverter.class)
			//			int max_cps;

			@Parameter(names="-internal", description="whether or not to identify internal change points", converter=BooleanConverter.class)
			boolean internal;

			@Parameter(names="-coverage", description="whether or not to calculate coverage of resulting models", converter=BooleanConverter.class, arity=1)
			boolean coverage=true;

			@Parameter(names="-min_fold", description="the minimum fold change between neighboring segments expressed as a ratio of low/high, acceptable values in the range [0.0-1.0]", converter=DoubleConverter.class)
			double min_fold = 0.5;

			@Parameter(names="-min_terminal", description="terminal segments are \"virtually\" extended by this amount before segmentation", converter=IntegerConverter.class)
			int min_terminal = 300;

			@Parameter(names="-jnct_alpha", description="The significance level for binomial test to accept splice junction", converter=DoubleConverter.class)
			double jnct_alpha = 0.05;

			@Parameter(names="-merge_segments", description="Optional:bed file of regions to merge across")
			String filled_gap_segments;

		}

		class MultiSegmentCommand{	

			@Parameter(names="-bam1", description="bam from sample 1")
			String bam1;

			@Parameter(names="-bam2", description="bam from sample 2")
			String bam2;

			@Parameter(names="-base1", description="base id from the assembly step for sample 1")
			String base1;

			@Parameter(names="-base2", description="base id from the assembly step for sample 2")
			String base2;

			@Parameter(names="-out", description="gtf to which output will be written")
			String out;

			@Parameter(names="-s", description="strandedness")
			String strandedness;

			@Parameter(names="-dir", description="The output directory for the assembly step")
			public String dir;
		}

		class ListCommand{	

			@Parameter(names="-pairfile", description="gtf from the pairwise analysis")
			String pairfile;

		}

		class EnumerateCommand{	
			@Parameter(names="-splicegraph", description="splice graph gtf from the assembly step")
			File splicegraph;
			
			@Parameter(names="-max_isoforms", description="loci with more than this number of isoforms will be skipped", converter=IntegerConverter.class)
			Integer max_paths;
			
			@Parameter(names="-base", description="splice isoforms will be written to [base].isoforms.gtf, skipped locus IDs of [base].skipped.txt")
			String base;

		}

		class HelpCommand{	

		}

		JCommander jc = new JCommander();
		AssembleCommand assemble = new AssembleCommand();
		SegmentCommand segment = new SegmentCommand();
		MultiSegmentCommand compare = new MultiSegmentCommand();
		ListCommand list = new ListCommand();
		EnumerateCommand enumerate = new EnumerateCommand();
		HelpCommand help = new HelpCommand();
		jc.addCommand("assemble", assemble);
		jc.addCommand("segment", segment);
		jc.addCommand("compare", compare);
		jc.addCommand("enumerate", enumerate);
			jc.addCommand("list", list);
		jc.addCommand("-h", help);
		jc.parse(args);

		if(jc.getParsedCommand()==null || jc.getParsedCommand().equals("-h")){
			jc.usage();
		}
		else if(jc.getParsedCommand().equals("assemble")){
			//pre-process
			{

				File out_dir = new File(assemble.dir);

				File bamFile = assemble.bam;

				File tmp_dir = FileUtils.getFile(out_dir,"tmp");
				if(!tmp_dir.exists())
					tmp_dir.mkdirs();

				if(!tmp_dir.exists())
					tmp_dir.mkdirs();

				File splice_junction_bed = FileUtils.getFile(tmp_dir,assemble.base+".sj.bed");
				File splice_count_gtf = FileUtils.getFile(tmp_dir,assemble.base+".sj.counted.gtf");
				File segment_bed = FileUtils.getFile(tmp_dir,assemble.base+".seg.bed");

				Strandedness strandedness = Strandedness.valueOf(assemble.strandedness);

				SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
				SAMFileReader sfr = new SAMFileReader(bamFile);
				System.out.println("finding junctions");
				BEDWriter sj_bw = new BEDWriter(IO.bufferedPrintstream(splice_junction_bed));
				FindSpliceJunctions.tabulateSpliceJunctions(sfr, sj_bw);
				sj_bw.close();

				System.out.println("counting junction supporting reads");
				GTFWriter counted_sj_gtf = new GTFWriter(IO.bufferedPrintstream(splice_count_gtf));
				FindSpliceJunctions.countJunctionSupportingReads(sfr, strandedness, splice_junction_bed, counted_sj_gtf);
				counted_sj_gtf.close();

				System.out.println("finding expressed segments");
				BEDWriter seg_bw = new BEDWriter(IO.bufferedPrintstream(segment_bed));
				SlidingWindow.identifyExpressed(sfr, strandedness, 1, assemble.threshold, 1, seg_bw);
				seg_bw.close();

				if(assemble.insert_size_quantile!=null){
					File matepair_bed = FileUtils.getFile(tmp_dir,assemble.base+".mate_pair.bed");
					File scaffolded_bed = FileUtils.getFile(tmp_dir,assemble.base+".scaffolded.bed");
					System.out.println("calculating insert size quantile");
					int insert_size = ClusterExpressedSegments.estimateMateInsertSize(sfr, strandedness, segment_bed, assemble.insert_size_quantile);
					System.out.printf("%.3fth-quantile insert size = %s bp\n", assemble.insert_size_quantile, insert_size);
					System.out.println("identifying proper mate pairs");
					BEDWriter matepair_bw = new BEDWriter(IO.bufferedPrintstream(matepair_bed));
					ClusterExpressedSegments.identifySpannableRegions(sfr, strandedness, matepair_bw, insert_size);
					matepair_bw.close();
					System.out.println("scaffolding segments spanned by mate pairs");
					BEDWriter scaffolded_bw = new BEDWriter(IO.bufferedPrintstream(scaffolded_bed));
					ClusterExpressedSegments.scaffoldSpannableRegions(segment_bed, matepair_bed, scaffolded_bw);
					scaffolded_bw.close();
					//					scaffolded_bed.renameTo(segment_bed);
				}

				//output models
				sfr.close();
			}
			// segmentation
			{
				File assembly_gtf = new File(Util.sprintf("%s/%s.gtf",assemble.dir, assemble.base));
				File unspliced = new File(Util.sprintf("%s/%s.unspliced.gtf",assemble.dir, assemble.base));

				if(assembly_gtf.getParentFile()!=null && !assembly_gtf.getParentFile().exists())
					assembly_gtf.getParentFile().mkdirs();

				File bamFile = assemble.bam;

				File out_dir = new File(assemble.dir);

				File tmp_dir = FileUtils.getFile(out_dir,"tmp");
				if(!tmp_dir.exists())
					tmp_dir.mkdirs();

				// files from the pre-process			
				File splice_junction_bed = FileUtils.getFile(tmp_dir,assemble.base+".sj.bed");
				File splice_count_gtf = FileUtils.getFile(tmp_dir,assemble.base+".sj.counted.gtf");
				File segment_bed = assemble.insert_size_quantile==null ? FileUtils.getFile(tmp_dir,assemble.base+".seg.bed") : FileUtils.getFile(tmp_dir,assemble.base+".scaffolded.bed");

				double jnct_alpha = assemble.jnct_alpha;


				// files to write
				File merged_segment_bed = FileUtils.getFile(tmp_dir,assemble.base+".seg.merged.bed");
				File inferred_strand_segment_bed = FileUtils.getFile(tmp_dir,assemble.base+".seg.merged.inferred_strand.bed");
				File spliced_exon_gtf = FileUtils.getFile(tmp_dir,assemble.base+".exon.gtf");
				File acc_jnct_gtf = FileUtils.getFile(tmp_dir,assemble.base+".sj.acc.gtf");
				File rej_jnct_gtf = FileUtils.getFile(tmp_dir,assemble.base+".sj.rej.gtf");
				File trimmed_exon_gtf = FileUtils.getFile(tmp_dir,assemble.base+".exon.trimmed.gtf");
				File intronic_exon_gtf = FileUtils.getFile(tmp_dir,assemble.base+".intronic.gtf");
				File changepoint_gtf = FileUtils.getFile(tmp_dir,assemble.base+".cp.gtf");
				File filtered_changepoint_gtf = FileUtils.getFile(tmp_dir,assemble.base+".cp.filtered.gtf");
				File changepoint_exon_gtf = FileUtils.getFile(tmp_dir,assemble.base+".cp_exon.gtf");

				Strandedness strandedness = Strandedness.valueOf(assemble.strandedness);

				SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
				SAMFileReader sfr = new SAMFileReader(bamFile);

				System.out.println("merging gapped segments");
				BEDWriter merged_segments = new BEDWriter(IO.bufferedPrintstream(merged_segment_bed));
				File merge_regions = assemble.filled_gap_segments == null? null : new File(assemble.filled_gap_segments);
				ClusterExpressedSegments.mergeSegments(segment_bed, merge_regions, merged_segments, assemble.merge_radius);
				merged_segments.close();

				System.out.println("filtering splice junctions");
				GTFWriter passed_jnct_gw = new GTFWriter(IO.bufferedPrintstream(acc_jnct_gtf));
				GTFWriter failed_jnct_gw = new GTFWriter(IO.bufferedPrintstream(rej_jnct_gtf));
				FindSpliceJunctions.filterCountedJunctions(splice_junction_bed, splice_count_gtf, jnct_alpha, passed_jnct_gw, failed_jnct_gw);
				passed_jnct_gw.close();
				failed_jnct_gw.close();

				if(strandedness==Strandedness.unstranded){
					System.out.println("inferring strand from spliced reads");
					BEDWriter inferred_strand_segment = new BEDWriter(IO.bufferedPrintstream(inferred_strand_segment_bed));
					ClusterExpressedSegments.inferSegmentStrand(acc_jnct_gtf, merged_segment_bed, inferred_strand_segment);
					inferred_strand_segment.close();
					merged_segment_bed = inferred_strand_segment_bed;
				}

				//cluster
				System.out.println("identifying spliced exons");
				//				double expressed_threshold = segment.threshold;
				GTFWriter spliced_exons = new GTFWriter(IO.bufferedPrintstream(spliced_exon_gtf));
				ClusterExpressedSegments.identifySplicedExons(acc_jnct_gtf, merged_segment_bed, spliced_exons);
				spliced_exons.close();

				if(strandedness==Strandedness.unstranded){
					System.out.println("trimming inferred exons that overlap on opposite strands");
					GTFWriter trimmed_exons = new GTFWriter(IO.bufferedPrintstream(trimmed_exon_gtf));
					ClusterExpressedSegments.trimInferredExons(spliced_exon_gtf, trimmed_exons);
					trimmed_exons.close();
					spliced_exon_gtf=trimmed_exon_gtf;
				}

				//				double expressed_threshold = segment.threshold;
				GTFWriter intronic_exon_writer = new GTFWriter(IO.bufferedPrintstream(intronic_exon_gtf));
				if(false){
					System.out.println("identifying intronic terminal exons");
					ClusterExpressedSegments.identifyIntronicExons(merged_segment_bed, acc_jnct_gtf, spliced_exon_gtf, sfr, strandedness, intronic_exon_writer, 10, .2);
				}
				intronic_exon_writer.close();

				// run SCM to identify changepoints in exons that are sufficiently long, use intronic and spliced exons as input.
				System.out.println("identifying change points");

				//negative binomial hyper parameters
				double alpha_0 = 1.0;
				double beta_0 = 1.0;
				int nb_r = assemble.nb_r;

				int r = assemble.segment_r; 
				double p = assemble.segment_p;

				int maxBins=20000;

				int binSize=assemble.w;
				int minCP=0;
				int minLength = binSize*2;
				GTFWriter changepoint_writer = new GTFWriter(IO.bufferedPrintstream(changepoint_gtf));
				boolean internal=assemble.internal;
				double min_fold = assemble.min_fold;
				int min_terminal = -assemble.min_terminal;
				IdentifyChangePoints.identifyNegativeBinomialChangePointsInLongSegments(sfr, spliced_exon_gtf, intronic_exon_gtf, changepoint_writer, minLength, maxBins, binSize, minCP, strandedness, alpha_0, beta_0, nb_r, r, p, internal, min_fold, min_terminal);
				changepoint_writer.close();

				System.out.println("filtering change points");
				int sj_radius = binSize*1;
				GTFWriter filtered_changepoint_writer = new GTFWriter(IO.bufferedPrintstream(filtered_changepoint_gtf));
				IdentifyChangePoints.filterChangePoints(sfr, spliced_exon_gtf, intronic_exon_gtf, acc_jnct_gtf, changepoint_gtf, filtered_changepoint_writer, sj_radius);
				filtered_changepoint_writer.close();

				System.out.println("enumerating change point exons");
				GTFWriter changepoint_exon_writer = new GTFWriter(IO.bufferedPrintstream(changepoint_exon_gtf));
				IdentifyChangePoints.enumerateChangepointExons(spliced_exon_gtf, intronic_exon_gtf, filtered_changepoint_gtf, changepoint_exon_writer);
				changepoint_exon_writer.close();

				//output models
				System.out.println("identifying connected components");
				GTFWriter assembly_writer = new GTFWriter(IO.bufferedPrintstream(assembly_gtf));
				ExonSpliceGraph.labelConnectedComponents(acc_jnct_gtf, spliced_exon_gtf, intronic_exon_gtf, changepoint_exon_gtf, assembly_writer);
				assembly_writer.close();

				System.out.println("identifying unspliced segments");
				GTFWriter remainder_writer = new GTFWriter(IO.bufferedPrintstream(unspliced));
				ClusterExpressedSegments.writeUnsplicedRemainder(merged_segment_bed, assembly_gtf, remainder_writer);
				remainder_writer.close();

				if(assemble.coverage){
					System.out.println("calculating coverage");
					File coverage_gtf = new File(Util.sprintf("%s/%s.coverage.gtf",assemble.dir, assemble.base));
					File coverage_unspliced = new File(Util.sprintf("%s/%s.unspliced.coverage.gtf",assemble.dir, assemble.base));

					GTFWriter coverage_gw = new GTFWriter(IO.bufferedPrintstream(coverage_gtf));
					IdentifyChangePoints.calculate_spliced_coverage(sfr, assembly_gtf, strandedness, coverage_gw);
					coverage_gw.close();
					GTFWriter coverage_unspliced_gw = new GTFWriter(IO.bufferedPrintstream(coverage_unspliced));
					IdentifyChangePoints.calculate_coverage(sfr, unspliced, strandedness, coverage_unspliced_gw);
					coverage_unspliced_gw.close();
					//IdentifyChangePoints.enumerateChangepointExons(spliced_exon_gtf, intronic_exon_gtf, filtered_changepoint_gtf, changepoint_exon_writer);	
				}

				//output models
				sfr.close();

			}
		}
		else if(jc.getParsedCommand().equals("segment")){	
			// segmentation
			{
				File assembly_gtf = new File(Util.sprintf("%s/%s.gtf",segment.dir, segment.base));
				File unspliced = new File(Util.sprintf("%s/%s.unspliced.gtf",segment.dir, segment.base));

				if(assembly_gtf.getParentFile()!=null && !assembly_gtf.getParentFile().exists())
					assembly_gtf.getParentFile().mkdirs();

				File bamFile = segment.bam;

				File out_dir = new File(segment.dir);

				File tmp_dir = FileUtils.getFile(out_dir,"tmp");
				if(!tmp_dir.exists())
					tmp_dir.mkdirs();

				// files from the pre-process			
				File splice_junction_bed = FileUtils.getFile(tmp_dir,segment.base+".sj.bed");
				File splice_count_gtf = FileUtils.getFile(tmp_dir,segment.base+".sj.counted.gtf");
				File segment_bed = !FileUtils.getFile(tmp_dir,segment.base+".scaffolded.bed").exists() ? FileUtils.getFile(tmp_dir,segment.base+".seg.bed") : FileUtils.getFile(tmp_dir,segment.base+".scaffolded.bed");

				double jnct_alpha = segment.jnct_alpha;


				// files to write
				File merged_segment_bed = FileUtils.getFile(tmp_dir,segment.base+".seg.merged.bed");
				File inferred_strand_segment_bed = FileUtils.getFile(tmp_dir,segment.base+".seg.merged.inferred_strand.bed");
				File spliced_exon_gtf = FileUtils.getFile(tmp_dir,segment.base+".exon.gtf");
				File acc_jnct_gtf = FileUtils.getFile(tmp_dir,segment.base+".sj.acc.gtf");
				File rej_jnct_gtf = FileUtils.getFile(tmp_dir,segment.base+".sj.rej.gtf");
				File trimmed_exon_gtf = FileUtils.getFile(tmp_dir,segment.base+".exon.trimmed.gtf");
				File intronic_exon_gtf = FileUtils.getFile(tmp_dir,segment.base+".intronic.gtf");
				File changepoint_gtf = FileUtils.getFile(tmp_dir,segment.base+".cp.gtf");
				File filtered_changepoint_gtf = FileUtils.getFile(tmp_dir,segment.base+".cp.filtered.gtf");
				File changepoint_exon_gtf = FileUtils.getFile(tmp_dir,segment.base+".cp_exon.gtf");

				Strandedness strandedness = Strandedness.valueOf(segment.strandedness);

				SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
				SAMFileReader sfr = new SAMFileReader(bamFile);

				System.out.println("merging gapped segments");
				BEDWriter merged_segments = new BEDWriter(IO.bufferedPrintstream(merged_segment_bed));
				File merge_regions = segment.filled_gap_segments == null? null : new File(segment.filled_gap_segments);
				ClusterExpressedSegments.mergeSegments(segment_bed, merge_regions, merged_segments, segment.merge_radius);
				merged_segments.close();

				System.out.println("filtering splice junctions");
				GTFWriter passed_jnct_gw = new GTFWriter(IO.bufferedPrintstream(acc_jnct_gtf));
				GTFWriter failed_jnct_gw = new GTFWriter(IO.bufferedPrintstream(rej_jnct_gtf));
				FindSpliceJunctions.filterCountedJunctions(splice_junction_bed, splice_count_gtf, jnct_alpha, passed_jnct_gw, failed_jnct_gw);
				passed_jnct_gw.close();
				failed_jnct_gw.close();

				if(strandedness==Strandedness.unstranded){
					System.out.println("inferring strand from spliced reads");
					BEDWriter inferred_strand_segment = new BEDWriter(IO.bufferedPrintstream(inferred_strand_segment_bed));
					ClusterExpressedSegments.inferSegmentStrand(acc_jnct_gtf, merged_segment_bed, inferred_strand_segment);
					inferred_strand_segment.close();
					merged_segment_bed = inferred_strand_segment_bed;
				}

				//cluster
				System.out.println("identifying spliced exons");
				//				double expressed_threshold = segment.threshold;
				GTFWriter spliced_exons = new GTFWriter(IO.bufferedPrintstream(spliced_exon_gtf));
				ClusterExpressedSegments.identifySplicedExons(acc_jnct_gtf, merged_segment_bed, spliced_exons);
				spliced_exons.close();

				if(strandedness==Strandedness.unstranded){
					System.out.println("trimming inferred exons that overlap on opposite strands");
					GTFWriter trimmed_exons = new GTFWriter(IO.bufferedPrintstream(trimmed_exon_gtf));
					ClusterExpressedSegments.trimInferredExons(spliced_exon_gtf, trimmed_exons);
					trimmed_exons.close();
					spliced_exon_gtf=trimmed_exon_gtf;
				}

				//				double expressed_threshold = segment.threshold;
				GTFWriter intronic_exon_writer = new GTFWriter(IO.bufferedPrintstream(intronic_exon_gtf));
				if(false){
					System.out.println("identifying intronic terminal exons");
					ClusterExpressedSegments.identifyIntronicExons(merged_segment_bed, acc_jnct_gtf, spliced_exon_gtf, sfr, strandedness, intronic_exon_writer, 10, .2);
				}
				intronic_exon_writer.close();

				// run SCM to identify changepoints in exons that are sufficiently long, use intronic and spliced exons as input.
				System.out.println("identifying change points");

				//negative binomial hyper parameters
				double alpha_0 = 1.0;
				double beta_0 = 1.0;
				int nb_r = segment.nb_r;

				int r = segment.segment_r; 
				double p = segment.segment_p;

				int maxBins=20000;

				int binSize=segment.w;
				int minCP=0;
				int minLength = binSize*2;
				GTFWriter changepoint_writer = new GTFWriter(IO.bufferedPrintstream(changepoint_gtf));
				boolean internal=segment.internal;
				double min_fold = segment.min_fold;
				int min_terminal = -segment.min_terminal;
				IdentifyChangePoints.identifyNegativeBinomialChangePointsInLongSegments(sfr, spliced_exon_gtf, intronic_exon_gtf, changepoint_writer, minLength, maxBins, binSize, minCP, strandedness, alpha_0, beta_0, nb_r, r, p, internal, min_fold, min_terminal);
				changepoint_writer.close();

				System.out.println("filtering change points");
				int sj_radius = binSize*1;
				GTFWriter filtered_changepoint_writer = new GTFWriter(IO.bufferedPrintstream(filtered_changepoint_gtf));
				IdentifyChangePoints.filterChangePoints(sfr, spliced_exon_gtf, intronic_exon_gtf, acc_jnct_gtf, changepoint_gtf, filtered_changepoint_writer, sj_radius);
				filtered_changepoint_writer.close();

				System.out.println("enumerating change point exons");
				GTFWriter changepoint_exon_writer = new GTFWriter(IO.bufferedPrintstream(changepoint_exon_gtf));
				IdentifyChangePoints.enumerateChangepointExons(spliced_exon_gtf, intronic_exon_gtf, filtered_changepoint_gtf, changepoint_exon_writer);
				changepoint_exon_writer.close();

				//output models
				System.out.println("identifying connected components");
				GTFWriter assembly_writer = new GTFWriter(IO.bufferedPrintstream(assembly_gtf));
				ExonSpliceGraph.labelConnectedComponents(acc_jnct_gtf, spliced_exon_gtf, intronic_exon_gtf, changepoint_exon_gtf, assembly_writer);
				assembly_writer.close();

				System.out.println("identifying unspliced segments");
				GTFWriter remainder_writer = new GTFWriter(IO.bufferedPrintstream(unspliced));
				ClusterExpressedSegments.writeUnsplicedRemainder(merged_segment_bed, assembly_gtf, remainder_writer);
				remainder_writer.close();

				if(segment.coverage){
					System.out.println("calculating coverage");
					File coverage_gtf = new File(Util.sprintf("%s/%s.coverage.gtf",segment.dir, segment.base));
					File coverage_unspliced = new File(Util.sprintf("%s/%s.unspliced.coverage.gtf",segment.dir, segment.base));

					GTFWriter coverage_gw = new GTFWriter(IO.bufferedPrintstream(coverage_gtf));
					IdentifyChangePoints.calculate_spliced_coverage(sfr, assembly_gtf, strandedness, coverage_gw);
					coverage_gw.close();
					GTFWriter coverage_unspliced_gw = new GTFWriter(IO.bufferedPrintstream(coverage_unspliced));
					IdentifyChangePoints.calculate_coverage(sfr, unspliced, strandedness, coverage_unspliced_gw);
					coverage_unspliced_gw.close();
					//IdentifyChangePoints.enumerateChangepointExons(spliced_exon_gtf, intronic_exon_gtf, filtered_changepoint_gtf, changepoint_exon_writer);	
				}

				//output models
				sfr.close();

			}
		}
		else if(jc.getParsedCommand().equals("compare")){
			String spliced_exon_gtf1 = compare.dir+"/tmp/"+compare.base1+".exon.gtf";
			String spliced_exon_gtf2 = compare.dir+"/tmp/"+compare.base2+".exon.gtf";			
			JointSegmentation.performJointSegmentation(compare.base1, compare.base2, spliced_exon_gtf1, spliced_exon_gtf2, compare.bam1, compare.bam2, compare.out, Strandedness.valueOf(compare.strandedness));	
		}
		else if(jc.getParsedCommand().equals("list")){
			JointSegmentation.identifyBypassedRegions(new File(list.pairfile));	
		}
		else if(jc.getParsedCommand().equals("enumerate")){
			GTFWriter gw = new GTFWriter(IO.bufferedPrintstream(Util.sprintf("%s.isoforms.gtf", enumerate.base)));
			PrintStream skipped = IO.bufferedPrintstream(Util.sprintf("%s.skipped.txt", enumerate.base));
			splicegraph.ExonSpliceGraph.iterateSpliceIsoforms(enumerate.splicegraph,gw,skipped,enumerate.max_paths);
			gw.close();
			skipped.close();
		}
	}
}
