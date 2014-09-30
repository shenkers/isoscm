package executable;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import javax.xml.bind.JAXBException;

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

public class IsoSCM {

	public static void main(String[] args) throws IOException, JAXBException {

		
//		args = new String[]{
//				"assemble",
//				"-bam", "/home/sol/lailab/sol/mel_yak_vir/total_rna/Y/H.bam",
////				"-bam2", "/home/sol/lailab/sol/mel_yak_vir/total_rna/Y/T.bam",
////				"-base1", "YH",
////				"-base2", "YT",
////				"-out_base", "/home/sol/Y.H.T",
////				"-s", "reverse_forward",
////				"-dir", "/home/sol/lailab/sol/mel_yak_vir/isoscm/Y/"
//				};
		
		 
		JCommander jc = new JCommander();
		jc.setProgramName("java -jar IsoSCM.jar");
		
		AssembleCommand assemble = new AssembleCommand();
		SegmentCommand segment = new SegmentCommand();
		CompareCommand compare = new CompareCommand();
		EnumerateCommand enumerate = new EnumerateCommand();
		HelpCommand help = new HelpCommand();
		
		jc.addCommand("assemble", assemble);
		jc.addCommand("segment", segment);
		jc.addCommand("compare", compare);
		jc.addCommand("enumerate", enumerate);
		jc.addCommand("-h", help);
		jc.parse(args);

		if(jc.getParsedCommand()==null || jc.getParsedCommand().equals("-h")){
			jc.usage();
		}
		else if(jc.getParsedCommand().equals("assemble")){
			File out_dir = assemble.dir;
			if(!out_dir.exists())
				out_dir.mkdirs();

			File configuration_file = FileUtils.getFile(out_dir, Util.sprintf("%s.assembly_parameters.xml", assemble.base));
			ConfigurationIO.writeConfiguration(assemble, configuration_file);
			
			//pre-process
			{

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
		else if(jc.getParsedCommand().equals("enumerate")){
			AssembleCommand assemblyConfiguration = ConfigurationIO.readAssemblyConfiguration(enumerate.assemblyXml);
			
			File out_dir = assemble.dir;
			File configuration_xml = FileUtils.getFile(out_dir, Util.sprintf("%s.enumerate_parameters.xml", assemble.base));
			ConfigurationIO.writeConfiguration(enumerate, configuration_xml);
			
			PrintStream gtfFile = IO.bufferedPrintstream(FileUtils.getFile(assemblyConfiguration.dir,Util.sprintf("%s.isoforms.gtf", assemblyConfiguration.base)));
			PrintStream skipped = IO.bufferedPrintstream(FileUtils.getFile(assemblyConfiguration.dir,Util.sprintf("%s.skipped_loci.txt", assemblyConfiguration.base)));
			File assembly_gtf = new File(Util.sprintf("%s/%s.gtf",assemblyConfiguration.dir, assemblyConfiguration.base));
			GTFWriter gw = new GTFWriter(gtfFile);
			splicegraph.ExonSpliceGraph.iterateSpliceIsoforms(assembly_gtf,gw,skipped,enumerate.max_paths);
			gw.close();
			skipped.close();
		}
		else if(jc.getParsedCommand().equals("compare")){
			if(!compare.dir.exists())
				compare.dir.mkdirs();
			
			File configuration_xml = FileUtils.getFile(compare.dir,Util.sprintf("%s.compare_parameters.xml", compare.base));
			ConfigurationIO.writeConfiguration(compare, configuration_xml);
			
			AssembleCommand assemblyConfiguration1 = ConfigurationIO.readAssemblyConfiguration(compare.assemblyXml1);
			AssembleCommand assemblyConfiguration2 = ConfigurationIO.readAssemblyConfiguration(compare.assemblyXml2);
			
			File spliced_exon_gtf1 = FileUtils.getFile(assemblyConfiguration1.dir,"tmp",assemblyConfiguration1.base+".exon.gtf");
			File spliced_exon_gtf2 = FileUtils.getFile(assemblyConfiguration1.dir,"tmp",assemblyConfiguration1.base+".exon.gtf");
			
			File table = FileUtils.getFile(compare.dir,Util.sprintf("%s.txt",compare.base));
			File gtf = FileUtils.getFile(compare.dir,Util.sprintf("%s.gtf",compare.base));
			
			JointSegmentation.performJointSegmentation(assemblyConfiguration1.base, assemblyConfiguration2.base, spliced_exon_gtf1, spliced_exon_gtf2, assemblyConfiguration1.bam, assemblyConfiguration2.bam, table, gtf, Strandedness.valueOf(assemblyConfiguration1.strandedness), Strandedness.valueOf(assemblyConfiguration2.strandedness));	
		}
	}
}
