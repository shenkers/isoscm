package tools;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import tools.ParseGTF.TranscriptIterator;
import util.IO;
import util.Util;

public class GTFTools {

	public static class GTFWriter{
		PrintStream out;

		public GTFWriter(PrintStream out) throws FileNotFoundException {
			this.out = out;	
		}

		public GTFWriter(String fileName) throws FileNotFoundException {
			out = IO.bufferedPrintstream(fileName);	
		}

		public void write(String chr, int start, int end, char strand){
			out.printf("%s\tsol\texon\t%d\t%d\t0.0\t%s\t.\t\n", chr,start,end,strand);
		}

		public void write(String type, String chr, int start, int end, char strand, String attributeString){
			out.printf("%s\tsol\t%s\t%d\t%d\t0.0\t%s\t.\t%s\n", chr,type,start,end,strand,attributeString);
		}

		public void close(){
			out.close();
		}
	}

	public static class AnnotationParser{
		public StrandedGenomicIntervalTree<Map<String,Object>> exons;
		public StrandedGenomicIntervalTree<Map<String,Object>> cdss;
		public StrandedGenomicIntervalTree<Map<String,Object>> introns;
		public StrandedGenomicIntervalTree<Map<String,Object>> utr5ps;
		public StrandedGenomicIntervalTree<Map<String,Object>> utr3ps;
		public StrandedGenomicIntervalTree<Map<String,Object>> genes;
		public StrandedGenomicIntervalTree<Map<String,Object>> transcripts;
		
		public AnnotationParser(Iterable<AnnotatedRegion> transcript_iterator) throws FileNotFoundException {
			
			exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
			cdss = new StrandedGenomicIntervalTree<Map<String,Object>>();

			for(AnnotatedRegion r : transcript_iterator){
				if(r.annotation.equals("exon")){
					exons.add(r.chr,r.start,r.end,r.strand, r.attributes);
				}
				if(r.annotation.equals("CDS")){
					cdss.add(r.chr,r.start,r.end,r.strand, r.attributes);
				}
			}

			transcripts = IntervalTools.aggregateIntervals(exons, Util.list("transcript_id","gene_id","gene_name"), "transcript_id");
			genes = IntervalTools.aggregateIntervals(transcripts, Util.list("gene_id","gene_name"), "gene_id");
			
		}

		public AnnotationParser(String gtfFile) throws FileNotFoundException {
			TranscriptIterator ti = new TranscriptIterator(gtfFile);
			exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
			cdss = new StrandedGenomicIntervalTree<Map<String,Object>>();
//			introns = new StrandedGenomicIntervalTree<Map<String,Object>>();
//			utr5ps = new StrandedGenomicIntervalTree<Map<String,Object>>();
//			utr3ps = new StrandedGenomicIntervalTree<Map<String,Object>>();

			for(AnnotatedRegion r : ti){
				if(r.annotation.equals("exon")){
					exons.add(r.chr,r.start,r.end,r.strand, r.attributes);
				}
				if(r.annotation.equals("CDS")){
					cdss.add(r.chr,r.start,r.end,r.strand, r.attributes);
				}
			}

			transcripts = IntervalTools.aggregateIntervals(exons, Util.list("transcript_id","gene_id","gene_name"), "transcript_id");
			genes = IntervalTools.aggregateIntervals(transcripts, Util.list("gene_id","gene_name"), "gene_id");
			
//			StrandedGenomicIntervalTree<Map<String,Object>> cds_aggregates = IntervalTools.aggregateIntervals(cdss, "transcript_id");
//
//			String attribute = "transcript_id";
//
//			for(AnnotatedRegion transcript : transcripts){
//				StrandedGenomicIntervalSet utrs = new StrandedGenomicIntervalSet();
//				StrandedGenomicIntervalSet introns = new StrandedGenomicIntervalSet();
//
//				utrs.add(transcript.chr, transcript.start, transcript.end, transcript.strand);
//				introns.add(transcript.chr, transcript.start, transcript.end, transcript.strand);
//
//				for(AnnotatedRegion exon : IntervalTools.SelectIntervals(exons, transcript.chr, transcript.start, transcript.end, transcript.strand, attribute, transcript.getAttribute(attribute))){
//					introns.remove(exon.chr, exon.start, exon.end, exon.strand);
//				}
//
//				for(AnnotatedRegion intron : introns){
//					utrs.remove(intron.chr, intron.start, intron.end, intron.strand);
//					Map<String, Object> attributes = new HashMap<String, Object>();
//					attributes.put(attribute, transcript.getAttribute(attribute));
//					this.introns.add(intron.chr, intron.start, intron.end, intron.strand, attributes);
//				}
//
//				AnnotatedRegion cds_aggregate = null;
//				for(AnnotatedRegion CDS : IntervalTools.SelectIntervals(cds_aggregates, transcript.chr, transcript.start, transcript.end, transcript.strand, attribute, transcript.getAttribute(attribute))){
//					cds_aggregate = CDS;
//					utrs.remove(CDS.chr, CDS.start, CDS.end, CDS.strand);
//				}
//
//
//				if(cds_aggregate!=null){
//					for(AnnotatedRegion utr : utrs){
//						Map<String, Object> attributes = new HashMap<String, Object>();
//						attributes.put(attribute, transcript.getAttribute(attribute));
//						if(transcript.isNegativeStrand()^(cds_aggregate.get5Prime()-utr.get5Prime()<0))
//							utr3ps.add(utr.chr, utr.start, utr.end, utr.strand, attributes);
//						else
//							utr5ps.add(utr.chr, utr.start, utr.end, utr.strand, attributes);	
//					}
//				}
//			}

		}
		
		public AnnotationParser(String gtfFile, boolean exon_cds_only) throws FileNotFoundException {
			this(gtfFile);
			
			if(!exon_cds_only) {
				introns = new StrandedGenomicIntervalTree<Map<String,Object>>();
				utr5ps = new StrandedGenomicIntervalTree<Map<String,Object>>();
				utr3ps = new StrandedGenomicIntervalTree<Map<String,Object>>();
				
				String transcript_attribute = "transcript_id";
				String gene_attribute = "gene_id";
				StrandedGenomicIntervalTree<Map<String,Object>> cds_aggregates = IntervalTools.aggregateIntervals(cdss, "transcript_id");
				
				for(AnnotatedRegion transcript : transcripts){
					StrandedGenomicIntervalSet utrs = new StrandedGenomicIntervalSet();
					StrandedGenomicIntervalSet introns = new StrandedGenomicIntervalSet();

					utrs.add(transcript.chr, transcript.start, transcript.end, transcript.strand);
					introns.add(transcript.chr, transcript.start, transcript.end, transcript.strand);

					for(AnnotatedRegion exon : IntervalTools.SelectIntervals(exons, transcript.chr, transcript.start, transcript.end, transcript.strand, transcript_attribute, transcript.getAttribute(transcript_attribute))){
						introns.remove(exon.chr, exon.start, exon.end, exon.strand);
					}

					for(AnnotatedRegion intron : introns){
						utrs.remove(intron.chr, intron.start, intron.end, intron.strand);
						Map<String, Object> attributes = new HashMap<String, Object>();
						attributes.put(transcript_attribute, transcript.getAttribute(transcript_attribute));
						attributes.put(gene_attribute, transcript.getAttribute(gene_attribute));
						this.introns.add(intron.chr, intron.start, intron.end, intron.strand, attributes);
					}

					AnnotatedRegion cds_aggregate = null;
					for(AnnotatedRegion CDS : IntervalTools.SelectIntervals(cds_aggregates, transcript.chr, transcript.start, transcript.end, transcript.strand, transcript_attribute, transcript.getAttribute(transcript_attribute))){
						cds_aggregate = CDS;
						utrs.remove(CDS.chr, CDS.start, CDS.end, CDS.strand);
					}


					if(cds_aggregate!=null){
						for(AnnotatedRegion utr : utrs){
							Map<String, Object> attributes = new HashMap<String, Object>();
							attributes.put(transcript_attribute, transcript.getAttribute(transcript_attribute));
							attributes.put(gene_attribute, transcript.getAttribute(gene_attribute));
								if(transcript.isNegativeStrand()^(cds_aggregate.get5Prime()-utr.get5Prime()<0))
								utr3ps.add(utr.chr, utr.start, utr.end, utr.strand, attributes);
							else
								utr5ps.add(utr.chr, utr.start, utr.end, utr.strand, attributes);	
						}
					}
				}
			}
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		AnnotationParser ap = new AnnotationParser("/home/sol/data/gene_models/Drosophila_melanogaster.BDGP5.69.gtf");
		for(AnnotatedRegion r : ap.utr3ps){
			System.out.println(r);
		}
	}
}
