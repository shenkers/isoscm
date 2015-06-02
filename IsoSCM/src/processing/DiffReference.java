package processing;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;

import tools.AnnotatedRegion;
import tools.IntervalTools;
import tools.StrandedGenomicIntervalTree;
import tools.GTFTools.AnnotationParser;
import tools.ParseGTF.TranscriptIterator;
import util.Util;

public class DiffReference {

	public static void labelExons(AnnotationParser annotation){
		
		for(AnnotatedRegion transcript : annotation.transcripts){
			List<AnnotatedRegion> exons = IntervalTools.SelectIntervals(annotation.exons, transcript.chr, transcript.start, transcript.end, transcript.strand, "transcript_id", transcript.getAttribute("transcript_id"));

			if(exons.size()>1){
				Collections.sort(exons,IntervalTools.AnnotatedRegionComparator(true, !transcript.isNegativeStrand()));
				List<AnnotatedRegion> internal = exons.subList(1, exons.size()-1);
				for(AnnotatedRegion r : internal){
					r.addAttribute("exon_type", "internal_exon");
				}
				exons.get(0).addAttribute("exon_type", "5p_exon");
				exons.get(exons.size()-1).addAttribute("exon_type", "3p_exon");
			}
			else{
				exons.get(0).addAttribute("exon_type", "unspliced");
			}
		}

		
	}


	public static void diff(String gtfFile, File compareGtf) throws FileNotFoundException {
		TranscriptIterator loci = new TranscriptIterator(compareGtf);
		AnnotationParser ap = new AnnotationParser(gtfFile);

		labelExons(ap);
		
		PrintStream out = System.out;
		
		out.printf("exon\tlocus_id\tmatch_5p_3p\tmatch_5p_3p_type\tmatch_5p\tmatch_5p_type\tmatch_3p\tmatch_3p_type\toverlaps\toverlaps_type\n");
		for(AnnotatedRegion r : loci){
			if(r.annotation.equals("locus")){
				String locus_id = (String) r.getAttribute("id");
				
				Set<String> matched5p = new HashSet<String>();
				Set<String> matched3p = new HashSet<String>();
				Set<String> matched5p3p = new HashSet<String>();
				Set<String> overlappingIds = new HashSet<String>();
				Set<String> overlappingTypes = new HashSet<String>();
				Set<String> matched5pTypes= new HashSet<String>();
				Set<String> matched3pTypes= new HashSet<String>();
				Set<String> matched5p3pTypes = new HashSet<String>();
				for(AnnotatedRegion e : ap.exons.overlappingRegions(r.chr, r.start, r.end, r.strand)){
					String transcriptId = (String) e.getAttribute("transcript_id");
					String exon_type = (String) e.getAttribute("exon_type");
					
					overlappingIds.add(transcriptId);
					overlappingTypes.add(exon_type);
					
					if(r.get5Prime()==e.get5Prime()){
						matched5p.add(transcriptId);
						matched5pTypes.add(exon_type);
					}
					if(r.get3Prime()==e.get3Prime()){
						matched3p.add(transcriptId);
						matched3pTypes.add(exon_type);
					}
					if(r.get5Prime()==e.get5Prime() && r.get3Prime()==e.get3Prime()){
						matched5p3p.add(transcriptId);
						matched5p3pTypes.add(exon_type);
					}
				}		
				
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", r, locus_id,
						StringUtils.join(matched5p3p,","),
						StringUtils.join(matched5p3pTypes,","),
						StringUtils.join(matched5p,","),
						StringUtils.join(matched5pTypes,","),
						StringUtils.join(matched3p,","),
						StringUtils.join(matched3pTypes,","),
						StringUtils.join(overlappingIds,","),
						StringUtils.join(overlappingTypes,",")
						);
			}
		}

	}
}
