package splicegraph;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Stack;

import tools.AnnotatedRegion;
import tools.GTFTools.GTFWriter;
import tools.IntervalTools;
import tools.ParseGTF.TranscriptIterator;
import tools.StrandedGenomicIntervalSet;
import tools.StrandedGenomicIntervalTree;
import util.Util;
import util.Util.ExtremeTracker;
import util.Util.MapCounter;

public class ExonSpliceGraph {

	public static class ConnectedComponentResult{
		public StrandedGenomicIntervalTree<Map<String, Object>> exons;
		public StrandedGenomicIntervalTree<Map<String, Object>> splice_junctions;

		public ConnectedComponentResult(StrandedGenomicIntervalTree<Map<String, Object>> exons, StrandedGenomicIntervalTree<Map<String, Object>> splice_junctions) {
			this.exons = exons;
			this.splice_junctions = splice_junctions;
		}
	}

	public static void labelConnectedComponents(File acc_jnct_gtf, File spliced_exon_gtf, File intronic_exon_gtf, File changepoint_exon_gtf, GTFWriter assembly_writer) throws FileNotFoundException{

		StrandedGenomicIntervalTree<Map<String,Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String,Object>> sj = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalSet exon_intervals = new StrandedGenomicIntervalSet();
		StrandedGenomicIntervalTree<Map<String,Object>> exon_interval_tree = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalSet clusters = new StrandedGenomicIntervalSet();

		TranscriptIterator splice_junctions = new TranscriptIterator(acc_jnct_gtf);
		for(AnnotatedRegion splice_junction : splice_junctions){
			if(!sj.contains(splice_junction.chr, splice_junction.start, splice_junction.end, splice_junction.strand)){
				sj.add(splice_junction.chr, splice_junction.start, splice_junction.end, splice_junction.strand, new HashMap<String, Object>());
				clusters.add(splice_junction.chr, splice_junction.start, splice_junction.end, splice_junction.strand);
			}
		}		


		for(File gtf : Util.list(spliced_exon_gtf, intronic_exon_gtf, changepoint_exon_gtf)){
			TranscriptIterator ti = new TranscriptIterator(gtf);
			for(AnnotatedRegion r : ti){
				//				if(!exons.contains(r.chr, r.start, r.end, r.strand)){
				exons.add(r.chr, r.start, r.end, r.strand, r.attributes);
				exon_intervals.add(r.chr, r.start, r.end, r.strand);
				clusters.add(r.chr, r.start, r.end, r.strand);
				//				}
			}
		}

		for(AnnotatedRegion exon_interval : exon_intervals){
			exon_interval_tree.add(exon_interval.chr, exon_interval.start, exon_interval.end, exon_interval.strand, new HashMap<String,Object>());
		}


		int componentCount = -1;
		for(AnnotatedRegion cluster : clusters){
			for(AnnotatedRegion exon_interval : exon_interval_tree.overlappingRegions(cluster.chr, cluster.start, cluster.end, cluster.strand)){
				Integer label = (Integer) exon_interval.getAttribute("label");
				if(label==null){
					componentCount++;
					exon_interval.addAttribute("label", componentCount);

					Queue<AnnotatedRegion> connected_intervals = new LinkedList<AnnotatedRegion>();
					connected_intervals.add(exon_interval);

					do{
						AnnotatedRegion current_interval = connected_intervals.poll();

						// find all junctions that connect two intervals
						for(AnnotatedRegion splice_junction : sj.overlappingRegions(current_interval.chr, current_interval.start-1, current_interval.end+1, current_interval.strand)){

							if((splice_junction.start >= current_interval.start-1 && splice_junction.start <= current_interval.end+1) || (splice_junction.end >= current_interval.start-1 && splice_junction.end <= current_interval.end+1))
								splice_junction.addAttribute("label", componentCount);


							// if the junction connects to another interval

							if(splice_junction.start >= current_interval.start-1 && splice_junction.end > current_interval.end+1){
								//								splice_junction.addAttribute("label", componentCount);

								for(AnnotatedRegion other_interval : exon_interval_tree.overlappingRegions(splice_junction.chr, splice_junction.end+1, splice_junction.end+1, splice_junction.strand)){
									Integer other_label = (Integer) other_interval.getAttribute("label");
									if(other_label==null){
										other_interval.addAttribute("label", componentCount);
										connected_intervals.add(other_interval);
									}
								}
							}

							if(splice_junction.end <= current_interval.end+1 && splice_junction.start < current_interval.start-1){
								//								splice_junction.addAttribute("label", componentCount);

								for(AnnotatedRegion other_interval : exon_interval_tree.overlappingRegions(splice_junction.chr, splice_junction.start-1, splice_junction.start-1, splice_junction.strand)){
									Integer other_label = (Integer) other_interval.getAttribute("label");
									if(other_label==null){
										other_interval.addAttribute("label", componentCount);
										connected_intervals.add(other_interval);
									}
								}
							}
						}
					}while(!connected_intervals.isEmpty());
				}
			}
		}

		for(AnnotatedRegion r : exon_interval_tree){
			for(AnnotatedRegion exon : exons.overlappingRegions(r.chr, r.start, r.end, r.strand)){
				exon.addAttribute("locus_id", Util.sprintf("locus.%08d", r.getAttribute("label")));
				assembly_writer.write("exon", exon.chr, exon.start, exon.end, exon.strand, AnnotatedRegion.GTFAttributeString(exon.attributes));
			}

		}
		for(AnnotatedRegion r : sj){
			r.addAttribute("locus_id", Util.sprintf("locus.%08d", r.getAttribute("label")));
			r.attributes.remove("label");
			r.addAttribute("color", "#EE0000");
			assembly_writer.write("splice_jnct", r.chr, r.start, r.end, r.strand, AnnotatedRegion.GTFAttributeString(r.attributes));
		}

		/*
		 * rather than build graph of connected components, do a traversal of the exons, labeling as you go
		 * and writing those to which a label is just added
		 */
		/* WRITE ALL THE EXONS IN THE LOCUS
		int componentCount = -1;
		for(AnnotatedRegion cluster : clusters){
			List<ConnectedComponentResult> connectedComponents = ExonSpliceGraph.findConnectedComponents(exon, sj, cluster);

			for(ConnectedComponentResult component : connectedComponents){
				componentCount++;
				for(AnnotatedRegion e : component.exons){
					Map<String,Object> attributes = new HashMap<String, Object>();
					attributes.put("locus_id", Util.sprintf("locus.%08d", componentCount));
					attributes.putAll(e.attributes);
					assembly.write("exon", e.chr, e.start, e.end, e.strand, AnnotatedRegion.GTFAttributeString(attributes));
				}
				for(AnnotatedRegion splice_junct : component.splice_junctions){
					Map<String,Object> attributes = new HashMap<String, Object>();
					attributes.put("locus_id", Util.sprintf("locus.%08d", componentCount));
					assembly.write("splice_junction", splice_junct.chr, splice_junct.start, splice_junct.end, splice_junct.strand, AnnotatedRegion.GTFAttributeString(attributes));
				}
			}
		}
		 */
	}

	public static void identifyInvariableExons(File assembly_gtf, GTFWriter isoform_gtf) throws FileNotFoundException{
		TranscriptIterator ti = new TranscriptIterator(assembly_gtf);

		StrandedGenomicIntervalTree<Map<String, Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalSet exonic_intervals = new StrandedGenomicIntervalSet();
		StrandedGenomicIntervalTree<Map<String, Object>> sj = new StrandedGenomicIntervalTree<Map<String,Object>>();

		for(AnnotatedRegion r : ti){
			if(r.annotation.equals("exon")){
				exonic_intervals.add(r.chr, r.start, r.end, r.strand);
				exons.add(r.chr, r.start, r.end, r.strand, r.attributes);
			}
			if(r.annotation.equals("splice_jnct")){
				sj.add(r.chr, r.start, r.end, r.strand, r.attributes);
			}
		}

		StrandedGenomicIntervalSet invariable_intervals = new StrandedGenomicIntervalSet();
		for(AnnotatedRegion r : sj){
			int n=0;
			for(AnnotatedRegion overlap : sj.overlappingRegions(r.chr, r.start, r.end, r.strand)){
				n++;
			}
			System.out.printf("n %d\n", n);
			if(n==1){
				ExtremeTracker<Integer> bounds = new ExtremeTracker<Integer>();
				int n5p=0, n3p=0;
				for(AnnotatedRegion exonic_interval : exonic_intervals.overlappingRegions(r.chr, r.start-1, r.start-1, r.strand)){
					bounds.put(exonic_interval.start);
					for(AnnotatedRegion exon : exons.overlappingRegions(exonic_interval.chr, exonic_interval.start, exonic_interval.end, exonic_interval.strand)){
						n5p++;
					}
				}
				for(AnnotatedRegion exonic_interval : exonic_intervals.overlappingRegions(r.chr, r.end+1, r.end+1, r.strand)){
					bounds.put(exonic_interval.end);
					for(AnnotatedRegion exon : exons.overlappingRegions(exonic_interval.chr, exonic_interval.start, exonic_interval.end, exonic_interval.strand)){
						n3p++;
					}
				}
				System.out.printf("%d,%d\n", n5p, n3p);
				if(n5p==1 && n3p==1){
					invariable_intervals.add(r.chr, bounds.getMin(), bounds.getMax(), r.strand);
				}
			}
		}

		int i=0; 

		for(AnnotatedRegion r : invariable_intervals){
			Map<String,Object> attributes = new HashMap<String, Object>();
			attributes.put("transcript_id", i);
			//			System.out.println(r);

			int n=0;
			for(AnnotatedRegion exon : exons.overlappingRegions(r.chr, r.start, r.end, r.strand)){
				n++;
				isoform_gtf.write("exon", exon.chr, exon.start, exon.end, exon.strand, AnnotatedRegion.GTFAttributeString(attributes));
			}
			i++;
			if(n>2){
				//				System.out.println(r);
				System.out.printf("%d\t%s\n", n, r);
			}
		}

	}    

	public static void outputSpliceGraph(File assembly_gtf, GTFWriter splicegraph_gtf) throws FileNotFoundException{
		TranscriptIterator ti = new TranscriptIterator(assembly_gtf);

		//		StrandedGenomicIntervalTree<Map<String, Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		//		StrandedGenomicIntervalSet exonic_intervals = new StrandedGenomicIntervalSet();
		//		StrandedGenomicIntervalTree<Map<String, Object>> sj = new StrandedGenomicIntervalTree<Map<String,Object>>();

		MapCounter<String> locus_counter = new MapCounter<String>();
		for(AnnotatedRegion r : ti){
//			System.out.println(r.toAttributeString());
			String locus_id = (String) r.getAttribute("locus_id");
			int i = locus_counter.get(locus_id);
			locus_counter.increment(locus_id);
			Map<String,Object> attributes = new HashMap<String, Object>();
			attributes.put("name", locus_id);
			attributes.put("transcript_id", locus_id + "." + i);
			
			if(r.annotation.equals("exon")){
				//				exonic_intervals.add(r.chr, r.start, r.end, r.strand);
				//				exons.add(r.chr, r.start, r.end, r.strand, r.attributes);
				if(r.getAttribute("type").equals("3p_exon")){
					int start = Math.min(IntervalTools.offsetPosition(r.get5Prime(), 1, r.isNegativeStrand(), false),r.get3Prime());
					int end = Math.max(IntervalTools.offsetPosition(r.get5Prime(), 1, r.isNegativeStrand(), false),r.get3Prime());
					splicegraph_gtf.write("exon", r.chr, start, end, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
				}
				if(r.getAttribute("type").equals("5p_exon")){
					int start = Math.min(IntervalTools.offsetPosition(r.get3Prime(), 1, r.isNegativeStrand(), true),r.get5Prime());
					int end = Math.max(IntervalTools.offsetPosition(r.get3Prime(), 1, r.isNegativeStrand(), true),r.get5Prime());
					splicegraph_gtf.write("exon", r.chr, start, end, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
				}
				if(r.getAttribute("type").equals("internal_exon")){
					int start = r.start+1;
					int end = r.end-1;
					splicegraph_gtf.write("exon", r.chr, start, end, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
				}
//				splicegraph_gtf.write("exon", r.chr, r.start, r.end, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
//				i++;
			}
			if(r.annotation.equals("splice_jnct")){
				splicegraph_gtf.write("exon", r.chr, r.start-1, r.start-1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
				splicegraph_gtf.write("exon", r.chr, r.end+1, r.end+1, r.strand, AnnotatedRegion.GTFAttributeString(attributes));
			}
		}

	}  

	public static void iterateSpliceIsoforms(File assembly_gtf, GTFWriter isoform_gtf) throws FileNotFoundException{
		TranscriptIterator ti = new TranscriptIterator(assembly_gtf);

		StrandedGenomicIntervalTree<Map<String, Object>> exon_5p = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> exon_3p = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> sj5p = new StrandedGenomicIntervalTree<Map<String,Object>>();

		MapCounter<String> isoform_count = new MapCounter<String>();

		for(AnnotatedRegion r : ti){
			if(r.annotation.equals("exon")){
				if(r.getAttribute("type").equals("5p_exon")){
					exon_5p.add(r.chr, r.start, r.end, r.strand, r.attributes);
				}
				else{
					exon_3p.add(r.chr, r.start, r.end, r.strand, r.attributes);
				}
			}
			if(r.annotation.equals("splice_jnct")){
				Map<String,Object> attributes = r.attributes;
				attributes.put("sj3p", r.get3Prime());
				sj5p.add(r.chr, r.get5Prime(), r.get5Prime(), r.strand, attributes);
			}
		}

		for(AnnotatedRegion r : exon_5p){	
			Stack<AnnotatedRegion> isoform_exons = new Stack<AnnotatedRegion>();
			iterateIsoforms(isoform_count, isoform_exons, r, exon_3p, sj5p, isoform_gtf);
		}
	}

	public static void iterateIsoforms(MapCounter<String> isoform_count, Stack<AnnotatedRegion> isoform_exons, AnnotatedRegion exon, StrandedGenomicIntervalTree<Map<String, Object>> exons, StrandedGenomicIntervalTree<Map<String, Object>> sj, GTFWriter gw){
		isoform_exons.push(exon);
		int i = exon.isNegativeStrand() ? -1 : 1;
		if(exon.getAttribute("type").equals("5p_exon")){

			// for each compatible exon
			// visit
			// increment id

			for(AnnotatedRegion jnct : sj.overlappingRegions(exon.chr, exon.get3Prime()+i, exon.get3Prime()+i, exon.strand)){
				int exon_5p = (Integer) jnct.getAttribute("sj3p") + i;
				for(AnnotatedRegion spliced_exon : exons.overlappingRegions(exon.chr, exon_5p, exon_5p, exon.strand)){
					if(spliced_exon.get5Prime() == exon_5p){
						iterateIsoforms(isoform_count, isoform_exons, spliced_exon, exons, sj, gw);
					}
				}
			}
		}
		else{			
			if(exon.getAttribute("type").equals("3p_exon")){
				String locus_id = (String) exon.getAttribute("locus_id");
				String isoform_id = locus_id + "." + isoform_count.get(locus_id);
				for(AnnotatedRegion r : isoform_exons){
					r.addAttribute("transcript_id", isoform_id);
					gw.write("exon", r.chr, r.start, r.end, r.strand, AnnotatedRegion.GTFAttributeString(r.attributes));
				}
				isoform_count.increment(locus_id);
			}
			else{
				for(AnnotatedRegion jnct : sj.overlappingRegions(exon.chr, exon.get3Prime()+i, exon.get3Prime()+i, exon.strand)){
					int exon_5p = (Integer) jnct.getAttribute("sj3p") + i;
					for(AnnotatedRegion spliced_exon : exons.overlappingRegions(exon.chr, exon_5p, exon_5p, exon.strand)){
						if(spliced_exon.get5Prime() == exon_5p){
							iterateIsoforms(isoform_count, isoform_exons, spliced_exon, exons, sj, gw);
						}
					}
				}
			}
		}
		isoform_exons.pop();
	}

	public static void main(String[] args) throws FileNotFoundException {
		{
			GTFWriter gw = new GTFWriter("/mnt/LaiLab/sol/mel_yak_vir/isoscm/M/MH.isoforms.gtf");
			iterateSpliceIsoforms(new File("/mnt/LaiLab/sol/mel_yak_vir/isoscm/M/MH.coverage.gtf"),gw);
			gw.close();
		}
		
		if(false)
		{
			GTFWriter gw = new GTFWriter("/home/sol/workspace/IsoSCM/tests/gtf/SRR594393.25.splice_graph.gtf");
			outputSpliceGraph(new File("/home/sol/workspace/IsoSCM/tests/gtf/SRR594393.25.gtf"), gw);
			gw.close();
		}
		System.exit(0);
		{
			GTFWriter gw = new GTFWriter("/home/sol/workspace/IsoSCM/tests/gtf/SRR594393.25.invariable.gtf");
			identifyInvariableExons(new File("/home/sol/workspace/IsoSCM/tests/gtf/SRR594393.25.gtf"), gw);
			gw.close();
		}
		System.exit(0);
		GTFWriter gw = new GTFWriter("/home/sol/workspace/IsoSCM/tests/gtf/SRR594393.25.isoforms.gtf");
		iterateSpliceIsoforms(new File("/home/sol/workspace/IsoSCM/tests/gtf/SRR594393.25.gtf"),gw);
		gw.close();
		System.exit(0);

		StrandedGenomicIntervalTree<Map<String, Object>> exons = new StrandedGenomicIntervalTree<Map<String,Object>>();
		StrandedGenomicIntervalTree<Map<String, Object>> splice_junctions = new StrandedGenomicIntervalTree<Map<String,Object>>();
		char strand = '+';
		String chr = "1";
		exons.add(chr, 1, 2, strand);
		exons.add(chr, 3, 4, strand);
		exons.add(chr, 5, 6, strand);
		exons.add(chr, 7, 8, strand);
		exons.add(chr, 9, 9, strand);
		splice_junctions.add(chr, 3, 6, strand);
		splice_junctions.add(chr, 5, 6, strand);
		splice_junctions.add(chr, 7, 8, strand);
//		UndirectedGraph<List<Integer>, DefaultEdge> g = createExonGraph(exons, splice_junctions, null);
//		System.out.println(labelConnectedComponents(g));
		System.exit(0);
		//		UndirectedGraph<List<Integer>, DefaultEdge> g = new SimpleGraph<List<Integer>, DefaultEdge>(DefaultEdge.class);
//		g.addVertex(Util.list(1,2));
//		g.addVertex(Util.list(1,2));
//		g.addVertex(Util.list(2,3));
//		g.addVertex(Util.list(3,3));
//		g.addVertex(Util.list(4,3));
//		g.addVertex(Util.list(4,4));
//		g.addVertex(Util.list(4,5));
//		g.addEdge(Util.list(1,2), Util.list(2,3));
//		g.addEdge(Util.list(1,2), Util.list(3,3));
//		g.addEdge(Util.list(3,3), Util.list(2,3));
//		g.addEdge(Util.list(4,3), Util.list(4,4));
//
//		DepthFirstIterator<List<Integer>, DefaultEdge> dfi = new DepthFirstIterator<List<Integer>, DefaultEdge>(g,Util.list(1,2));
//		while(dfi.hasNext()){
//			System.out.println(dfi.next());
//		}
//		System.out.println();
//		for(List<Integer> v : g.vertexSet()){
//			System.out.println(v);
//		}
	}
}
