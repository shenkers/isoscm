package tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import util.Util.ExtremeTracker;
import filter.ComposableFilter;
import filter.Counter;

public class BAMTools {

	public static class BAMWriter{
		SAMFileWriter sfw;
		SAMFileHeader sfh;

		/**
		 * Assumes the file is not sorted, and the desired sort order is coordinate
		 * @param file
		 * @param sfh
		 * @throws FileNotFoundException
		 */
		public BAMWriter(File file, SAMFileHeader sfh) throws FileNotFoundException {
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
			SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			this.sfh = sfh;
			sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);
			sfw = sfwf.makeBAMWriter(sfh, false, file);
		}
		
		/**
		 * Assumes the file is not sorted, and the desired sort order is coordinate
		 * @param file
		 * @param sfh
		 * @param presorted
		 * @throws FileNotFoundException
		 */
		public BAMWriter(File file, SAMFileHeader sfh, boolean presorted, boolean createIndex) throws FileNotFoundException {
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(createIndex);
			SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			this.sfh = sfh;
			sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);
			sfw = sfwf.makeBAMWriter(sfh, presorted, file);
		}

		/**
		 * Assumes the file is not sorted, and the desired sort order is coordinate
		 * @param file
		 * @param sfh
		 * @param presorted
		 * @throws FileNotFoundException
		 */
		public BAMWriter(File file, SAMFileHeader sfh, boolean presorted) throws FileNotFoundException {
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
			SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			this.sfh = sfh;
			sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);
			sfw = sfwf.makeBAMWriter(sfh, presorted, file);
		}

		/**
		 * @param file
		 * @param sfh
		 * @param presorted
		 * @throws FileNotFoundException
		 */
		public BAMWriter(File file, SAMFileHeader sfh, SAMFileHeader.SortOrder sortOrder, boolean presorted) throws FileNotFoundException {
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
			SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			this.sfh = sfh;
			sfh.setSortOrder(sortOrder);
			sfw = sfwf.makeBAMWriter(sfh, presorted, file);
		}

		/**
		 * 
		 * @param file
		 * @param sfh
		 * @param sortOrder
		 * @param presorted
		 * @throws FileNotFoundException
		 */
		public BAMWriter(File file, SAMFileHeader sfh, SAMFileHeader.SortOrder sortOrder, boolean presorted, boolean isBam) throws FileNotFoundException {
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
			SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			this.sfh = sfh;
			sfh.setSortOrder(sortOrder);
			if(isBam)
				sfw = sfwf.makeBAMWriter(sfh, presorted, file);
			else
				sfw = sfwf.makeSAMWriter(sfh, presorted, file);
		}

		public void write(SAMRecord sr){
			sfw.addAlignment(sr);
		}

		public void write(String name, String chr, int start, int end, boolean isNegativeStrand){
			SAMRecord sr = new SAMRecord(sfh);
			sr.setReferenceName(chr);
			sr.setMappingQuality(255);
			sr.setReadNegativeStrandFlag(isNegativeStrand);
			sr.setAlignmentStart(start);
			sr.setCigarString((end-start+1)+"M");
			sr.setReadName(name);
			sfw.addAlignment(sr);
		}

		public void write(String name, String chr, int start, String cigar, boolean isNegativeStrand){
			SAMRecord sr = new SAMRecord(sfh);
			sr.setReferenceName(chr);
			sr.setMappingQuality(255);
			sr.setReadNegativeStrandFlag(isNegativeStrand);
			sr.setAlignmentStart(start);
			sr.setCigarString(cigar);
			sr.setReadName(name);
			sfw.addAlignment(sr);
		}

		public void write(String name, String chr, int start, String cigar, Map<String,Object> attributes, boolean isNegativeStrand){
			SAMRecord sr = new SAMRecord(sfh);
			sr.setReferenceName(chr);
			sr.setMappingQuality(255);
			sr.setReadNegativeStrandFlag(isNegativeStrand);
			sr.setAlignmentStart(start);
			sr.setCigarString(cigar);
			sr.setReadName(name);
			for(String tag : attributes.keySet()){
				sr.setAttribute(tag, attributes.get(tag));
			}
			sfw.addAlignment(sr);
		}

		public void write(String name, String chr, int start, String cigar, Map<String,Object> attributes, boolean isNegativeStrand, String sequence, String quality){
			SAMRecord sr = new SAMRecord(sfh);
			sr.setReferenceName(chr);
			sr.setMappingQuality(255);
			sr.setReadNegativeStrandFlag(isNegativeStrand);
			sr.setAlignmentStart(start);
			sr.setCigarString(cigar);
			sr.setReadName(name);
			sr.setReadString(sequence);
			sr.setBaseQualityString(quality);
			for(String tag : attributes.keySet()){
				sr.setAttribute(tag, attributes.get(tag));
			}
			sfw.addAlignment(sr);
		}

		public void close(){
			sfw.close();
		}
	}

	public static class CumulativeReadsInterval{
		boolean contained;
		int start,end;
		int l;
		int[] cumulativeReads;

		public CumulativeReadsInterval(SamReader sfr, String chr, int start, int end, boolean contained){
			cumulativeReads = new int[end-start+1];
			l=cumulativeReads.length;
			this.contained=contained;
			this.start=start;
			this.end=end;
			SAMRecordIterator sri = sfr.query(chr, start, end, contained);
			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				//				System.out.println(Util.list(sr.getAlignmentStart(),sr.getAlignmentEnd()));
				//				
				for(int i=start; i<sr.getAlignmentEnd()+1 && i-start < l; i++){
					cumulativeReads[i-start]++;
				}

			}
			sri.close();
		}

		public int nAlignedReads(int start, int end){
			int n1=0,n2=0;

			n2 = cumulativeReads[end-this.start];
			n1 = cumulativeReads[start-this.start];

			//				return cumulativeReads[end-this.start]-cumulativeReads[start-this.start];

			return n1-n2;

		}
	}

	/**
	 * 
	 * @param sfr
	 * @return the number of aligned reads in the given archive
	 */
	public static int totalAlignedReads(SamReader sfr){
		if(!sfr.hasIndex()){
			throw new IllegalArgumentException("bam file must be indexed");
		}
		int aligned = 0;

		SAMFileHeader sfh = sfr.getFileHeader();
		int i=0;
		SAMSequenceRecord ssr = sfh.getSequence(i);
		while(ssr!=null){
			aligned += sfr.indexing().getIndex().getMetaData(i).getAlignedRecordCount();
			i++;
			ssr = sfh.getSequence(i);	
		}

		return aligned;
	}

	/**
	 * 
	 * @param sfr
	 * @return the number of aligned reads in the given archive
	 */
	public static int totalUnmappedReads(SamReader sfr){
		if(!sfr.hasIndex()){
			throw new IllegalArgumentException("bam file must be indexed");
		}
		int aligned = 0;

		SAMRecordIterator sri = sfr.queryUnmapped();
		while(sri.hasNext()){
			aligned++;
			sri.next();
		}
		sri.close();

		return aligned;
	}

	public static double FPKM(SamReader sfr, AnnotatedRegion ar, boolean contained){
		return FPKM(sfr, ar.chr, ar.start, ar.end, contained,true);
	}

	public static double FPKM(SamReader sfr, String chr, int start, int end, boolean contained){
		double nTotal = totalAlignedReads(sfr);
		double nAligned = nAlignedReads(sfr, chr, start, end, contained);
		return (nAligned*1000000.0*1000.0)/(nTotal*(end-start+1)*1.0);
	}

	public static double FPKM(SamReader sfr, String chr, int start, int end, boolean contained, boolean consuming, boolean uniquelyMapping){
		double nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		if(consuming)
			nAligned = nConsumingReads(sfr, chr, start, end, contained, uniquelyMapping);
		else
			nAligned = nAlignedReads(sfr, chr, start, end, contained);
		return (nAligned*1000000.0*1000.0)/(nTotal*(end-start+1)*1.0);
	}

	public static double FPKM(SamReader sfr, String chr, int start, int end, boolean contained, boolean consuming){
		double nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		if(consuming)
			nAligned = nConsumingReads(sfr, chr, start, end, contained);
		else
			nAligned = nAlignedReads(sfr, chr, start, end, contained);
		return (nAligned*1000000.0*1000.0)/(nTotal*(end-start+1)*1.0);
	}

	public static boolean expressedAboveThreshold(SamReader sfr, String chr, int start, int end, boolean contained, double threshold){
		double nTotal = totalAlignedReads(sfr);
		int minReads = (int) ((threshold * nTotal * (end-start+1))/(1000000*1000));
		return moreThanNAlignedReads(sfr, chr, start, end, contained, minReads);
	}

	public static double FPKM(SamReader sfr, List<AnnotatedRegion> exons, boolean contained) {
		int nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		int size=0;
		for(AnnotatedRegion exon : exons){
			nAligned += nAlignedReads(sfr, exon, contained);
			size += exon.size;
		}
		return (nAligned*1000000.0*1000)/(nTotal*1.0*size);
	}

	public static double FPKM(SamReader sfr, List<AnnotatedRegion> exons, boolean contained, boolean isNegativeStrand) {
		int nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		int size=0;
		for(AnnotatedRegion exon : exons){
			nAligned += nAlignedReads(sfr, exon.chr, exon.start, exon.end, contained, isNegativeStrand);
			size += exon.size;
		}
		return (nAligned*1000000.0*1000)/(nTotal*1.0*size);
	}

	public static double FPKM(SamReader sfr, List<AnnotatedRegion> exons, boolean contained, Strandedness strandedness, boolean isNegativeStrand) {
		int nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		int size=0;
		for(AnnotatedRegion exon : exons){
			nAligned += nConsumingReads(sfr, exon.chr, exon.start, exon.end, contained, strandedness, isNegativeStrand);
			size += exon.size;
		}
		return (nAligned*1000000.0*1000)/(nTotal*1.0*size);
	}

	public static double FPKM(SamReader sfr, AnnotatedRegion exon, Strandedness strandedness, boolean contained) {
		int nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		int size=0;

		nAligned += nConsumingReads(sfr, exon.chr, exon.start, exon.end, contained, strandedness, exon.isNegativeStrand());
		size += exon.size;

		return (nAligned*1000000.0*1000)/(nTotal*1.0*size);
	}

	public static double FPKM(SamReader sfr, String chr, int start, int end, boolean isNegativeStrand, boolean contained, Strandedness strandedness) {
		int nTotal = totalAlignedReads(sfr);
		int nAligned = 0;
		int size=end-start+1;

		nAligned += nConsumingReads(sfr, chr, start, end, contained, strandedness, isNegativeStrand);


		return (nAligned*1000000.0*1000)/(nTotal*1.0*size);
	}

	public static StrandedGenomicIntervalTree<Map<String,Object>> spliceJunctions(SAMRecord sr){
		StrandedGenomicIntervalTree<Map<String,Object>> jncts = new StrandedGenomicIntervalTree<Map<String,Object>>();
		if(sr.getCigarString().indexOf('N') == -1)
			return jncts;

		String chr = sr.getReferenceName();
		Cigar cigar = sr.getCigar();
		int alignmentPosition = sr.getAlignmentStart();

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			if(cigarElement.getOperator().consumesReferenceBases()){
				boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
				if(isIntron){
					AnnotatedRegion sj = new AnnotatedRegion("", chr, alignmentPosition,alignmentPosition+cigarElement.getLength()-1, sr.getCharacterAttribute("XS"));
					jncts.add(sj);
				}
				alignmentPosition += cigarElement.getLength();
			}			
		}

		return jncts;
	}
	
	public static StrandedGenomicIntervalTree<Map<String,Object>> consumingSegments(SAMRecord sr, Strandedness strandedness){
		StrandedGenomicIntervalTree<Map<String,Object>> segments = new StrandedGenomicIntervalTree<Map<String,Object>>();
		
		if(!sr.getReadUnmappedFlag()){
		String chr = sr.getReferenceName();
		Cigar cigar = sr.getCigar();
		int alignmentPosition = sr.getAlignmentStart();
		char strand = strand(sr, strandedness);

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			if(cigarElement.getOperator().consumesReferenceBases()){
				boolean consumesReference = cigarElement.getOperator().consumesReadBases();
				if(consumesReference){
					AnnotatedRegion sj = new AnnotatedRegion("", chr, alignmentPosition,alignmentPosition+cigarElement.getLength()-1, strand);
					segments.add(sj);
				}
				alignmentPosition += cigarElement.getLength();
			}			
		}
		}

		return segments;
	}

	public static GenomicIntervalTree<Map<String,Object>> deletions(SAMRecord sr){
		GenomicIntervalTree<Map<String,Object>> dels = new GenomicIntervalTree<Map<String,Object>>();
		if(sr.getCigarString().indexOf('D') == -1)
			return dels;

		String chr = sr.getReferenceName();
		Cigar cigar = sr.getCigar();
		int alignmentPosition = sr.getAlignmentStart();

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			if(cigarElement.getOperator().consumesReferenceBases()){
				boolean isDeletion = cigarElement.getOperator().equals(CigarOperator.D);
				if(isDeletion){
					dels.add(chr,alignmentPosition,alignmentPosition+cigarElement.getLength()-1);
				}
				alignmentPosition += cigarElement.getLength();
			}			
		}

		return dels;
	}

	public static int nAlignedSplicedReads(SamReader sfr, String chr, int start, int end, boolean contained){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			SAMRecord sr = sri.next();
			if(sr.getCigarString().indexOf('N') == -1)
				continue;

			Cigar cigar = sr.getCigar();
			int alignmentPosition = sr.getAlignmentStart();

			spliced:{
				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean isIntron = cigarElement.getOperator().equals(CigarOperator.N);
						if(isIntron && alignmentPosition >= start && alignmentPosition <= end){
							n++;
							break spliced;
						}
						for(int j=0; j<cigarElement.getLength(); j++){
							alignmentPosition++;
						}
						if(isIntron && alignmentPosition >= start && alignmentPosition <= end){
							n++;
							break spliced;
						}
					}			
				}
			}
		}
		sri.close();
		return n;
	}

	public static int nConsumingReads(SamReader sfr, String chr, int start, int end, boolean contained){
		return nConsumingReads(sfr, chr, start, end, contained, false);
	}

	public static int nConsumingReads(SamReader sfr, String chr, int start, int end, boolean contained, boolean unique){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			SAMRecord sr = sri.next();

			if(sr.getDuplicateReadFlag())
				continue;

			if(unique){
				if((Integer) sr.getAttribute("NH") > 1)
					continue;
			}

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

			if(
					//					sr.getMappingQuality()>=10 && 
					consumesOverlappingBases)
				n++;
		}

		sri.close();
		return n;
	}

	public static int nConsumingReads(SamReader sfr, String chr, int start, int end, boolean contained, boolean unique, boolean isNegativeStrand){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			SAMRecord sr = sri.next();

			if(sr.getReadNegativeStrandFlag()^isNegativeStrand)
				continue;

			if(sr.getDuplicateReadFlag())
				continue;

			if(unique){
				if((Integer) sr.getAttribute("NH") > 1)
					continue;
			}

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

			if(
					//					sr.getMappingQuality()>=10 && 
					consumesOverlappingBases)
				n++;
		}

		sri.close();
		return n;
	}

	public static int nConsumingReads(SamReader sfr, String chr, int start, int end, boolean contained, Strandedness strandedness, boolean isNegativeStrand, int threshold){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			SAMRecord sr = sri.next();
			if(strandedness.equals(Strandedness.unstranded)){
				if(consumesOverlappingBases(sr, start, end)){
					n++;
				}
			}
			else{
				char strand = strand(sr, strandedness);
				if((strand=='+'||strand=='-')&&( (strand=='-') == isNegativeStrand)){
					if(consumesOverlappingBases(sr, start, end)){
						n++;
					}
				}
			}
			if(n==threshold){
				break;
			}
		}
		sri.close();
		return n;
	}

	public static int[] nConsumingReads(SamReader sfr, String chr, int start, int end, boolean contained, Strandedness strandedness, int threshold, int max_hits){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int[] n = null;
		if(strandedness.equals(Strandedness.unstranded)){
			n = new int[1];
		}
		else{
			n = new int[2];
		}
		while(sri.hasNext()){

			SAMRecord sr = sri.next();
			int NH = sr.getIntegerAttribute("NH");
			if(NH <= max_hits){
				if(strandedness.equals(Strandedness.unstranded)){
					if(consumesOverlappingBases(sr, start, end)){
						n[0]++;
					}
				}
				else{
					char strand = strand(sr, strandedness);
					if((strand=='+'||strand=='-')){
						if(consumesOverlappingBases(sr, start, end)){
							if(strand=='-')
								n[1]++;
							else
								n[0]++;
						}
					}
				}
				if(strandedness.equals(Strandedness.unstranded)){
					if(n[0]>=threshold){
						break;
					}
				}
				else{
					if(n[0]>=threshold && n[1]>=threshold){
						break;
					}
				}
			}
		}
		sri.close();
		return n;
	}

	public static int[] nConsumingReads(SamReader sfr, String chr, int start, int end, boolean contained, Strandedness strandedness, int max_hits){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int[] n = null;
		if(strandedness.equals(Strandedness.unstranded)){
			n = new int[1];
		}
		else{
			n = new int[2];
		}
		while(sri.hasNext()){

			SAMRecord sr = sri.next();
			int NH = sr.getIntegerAttribute("NH");
			if(NH <= max_hits){
				if(strandedness.equals(Strandedness.unstranded)){
					if(consumesOverlappingBases(sr, start, end)){
						n[0]++;
					}
				}
				else{
					char strand = strand(sr, strandedness);
					if((strand=='+'||strand=='-')){
						if(consumesOverlappingBases(sr, start, end)){
							if(strand=='-')
								n[1]++;
							else
								n[0]++;
						}
					}
				}
			}
		}
		sri.close();
		return n;
	}

	public static int nConsumingReads(SamReader SAMReader, String chromosome, int start, int end, boolean contained, Strandedness strandedness, boolean isNegativeStrand){
		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

		int n=0;

		if(strandedness == Strandedness.unstranded){

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				if(consumesOverlappingBases(sr, start, end))
					n++;
			}

			sri.close();

			return n;
		}
		else if(strandedness == Strandedness.reverse_forward){

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();
				if(sr.getReadPairedFlag()){	
					if((sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag()) ^ isNegativeStrand){
						if(consumesOverlappingBases(sr, start, end))
							n++;
					}
				}
			}

			sri.close();

			return n;
		}
		else if(strandedness == Strandedness.forward){
			while(sri.hasNext() ){
				SAMRecord sr = sri.next();
				if(sr.getReadNegativeStrandFlag() == isNegativeStrand){
					if(consumesOverlappingBases(sr, start, end))
						n++;
				}
			}


			sri.close();

			return n;
		}
		else if(strandedness == Strandedness.reverse){
			while(sri.hasNext() ){
				SAMRecord sr = sri.next();
				if(sr.getReadNegativeStrandFlag() != isNegativeStrand){
					if(consumesOverlappingBases(sr, start, end))
						n++;
				}
			}


			sri.close();

			return n;
		}
		else{
			sri.close();
			throw new RuntimeException("Unimplemented");
		}
	}

	public static boolean allElementsConsumeReference(SAMRecord sr){
		Cigar cigar = sr.getCigar();

		boolean allConsumeReference = true;

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			allConsumeReference &= cigarElement.getOperator().consumesReferenceBases();
			if(!allConsumeReference)
				break;
		}

		return allConsumeReference;
	}
	
	public static boolean allElementsConsumeRead(SAMRecord sr){
		Cigar cigar = sr.getCigar();

		boolean allConsumeRead = true;

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			allConsumeRead &= cigarElement.getOperator().consumesReadBases();
			if(!allConsumeRead)
				break;
		}

		return allConsumeRead;
	}

	public static boolean consumesOverlappingBases(SAMRecord sr, int start, int end){
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

		return consumesOverlappingBases;
	}

	public static int[] terminalConsumedIndices(SAMRecord sr){
		ExtremeTracker<Integer> indices = new ExtremeTracker<Integer>();
		int alignmentPosition = sr.getAlignmentStart();
		int readPosition = 0;
		Cigar cigar = sr.getCigar();

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
			boolean consumesReferenceBases = cigarElement.getOperator().consumesReferenceBases();
			if(consumesReferenceBases && consumesReadBases){
				indices.put(readPosition);
				indices.put(readPosition+cigarElement.getLength()-1);
			}
			if(consumesReferenceBases)
				alignmentPosition += cigarElement.getLength();
			if(consumesReadBases)
				readPosition += cigarElement.getLength();

		}			



		return new int[]{indices.getMin(),indices.getMax()};
	}

	class FilteredReadCounter {

		ComposableFilter<SAMRecord> head;
		Counter<SAMRecord> counter;

		public FilteredReadCounter(List<ComposableFilter<SAMRecord>> f) {
			ComposableFilter<SAMRecord> s = f.get(0);
			head = s;
			for(int i=1; i<f.size(); i++){
				ComposableFilter<SAMRecord> n = f.get(1);
				s.append(n);
				s = n;
			}
			counter = new Counter<SAMRecord>();
			s.append(counter);
		}

		public int countFilteredReads(SAMRecordIterator sri){
			counter.reset();
			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				head.add(sr);
			}
			return counter.n;
		}
	}

	public static int nTerminatingReads(SamReader SAMReader, String chromosome, int pos, boolean end5p, Strandedness strandedness, boolean isNegativeStrand){
		SAMRecordIterator sri = SAMReader.query(chromosome, pos, pos, false);

		int n=0;

		if(strandedness == Strandedness.forward){

			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				if(sr.getReadNegativeStrandFlag() == isNegativeStrand){
					int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentStart() : sr.getAlignmentEnd();
					if(terminus == pos){
						n++;
					}
				}
			}

			sri.close();

			return n;
		}
		//		else if(strandedness == Strandedness.unstranded){
		//
		//			while(sri.hasNext() ){
		//				SAMRecord sr = sri.next();
		//
		//				if(consumesOverlappingBases(sr, start, end))
		//					n++;
		//			}
		//
		//			sri.close();
		//
		//			return n;
		//		}
		//		else if(strandedness == Strandedness.reverse_forward){
		//
		//			while(sri.hasNext() ){
		//				SAMRecord sr = sri.next();
		//				if(sr.getReadPairedFlag()){	
		//					if((sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag()) ^ isNegativeStrand){
		//						if(consumesOverlappingBases(sr, start, end))
		//							n++;
		//					}
		//				}
		//			}
		//
		//			sri.close();
		//
		//			return n;
		//		}
		//		else if(strandedness == Strandedness.reverse){
		//			while(sri.hasNext() ){
		//				SAMRecord sr = sri.next();
		//				if(sr.getReadNegativeStrandFlag() != isNegativeStrand){
		//					if(consumesOverlappingBases(sr, start, end))
		//						n++;
		//				}
		//			}
		//
		//
		//			sri.close();
		//
		//			return n;
		//		}
		else{
			sri.close();
			throw new RuntimeException("Unimplemented");
		}
	}

	/**
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param end
	 * @param end5p - choose the end of the read based on the transcript coordinates
	 * @param strandedness
	 * @param isNegativeStrand
	 * @return
	 */
	public static int nTerminatingReads(SamReader SAMReader, String chromosome, int start, int end, boolean end5p, Strandedness strandedness, boolean isNegativeStrand){
		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, false);

		int n=0;

		if(strandedness == Strandedness.forward){

			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				if(sr.getReadNegativeStrandFlag() == isNegativeStrand && sr.getIntegerAttribute("NH")==1){
					int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentStart() : sr.getAlignmentEnd();
					if(terminus >= start && terminus <= end){
						n++;
					}
				}
			}

			sri.close();

			return n;
		}
		if(strandedness == Strandedness.reverse){

			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				if(sr.getReadNegativeStrandFlag() != isNegativeStrand){
					int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentEnd() : sr.getAlignmentStart();
					if(terminus >= start && terminus <= end){
						n++;
					}
				}
			}

			sri.close();

			return n;
		}
		else{
			sri.close();
			throw new RuntimeException("Unimplemented");
		}
	}

	/**
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param end
	 * @param end5p - choose the end of the read based on the transcript coordinates
	 * @param strandedness
	 * @param isNegativeStrand
	 * @return
	 */
	public static int nTerminatingReads(SamReader SAMReader, String chromosome, int start, int end, boolean end5p, Strandedness strandedness, boolean isNegativeStrand, int max_hits){
		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, false);

		int n=0;

		if(strandedness == Strandedness.forward){

			while(sri.hasNext()){
				SAMRecord sr = sri.next();

				int NH = sr.getIntegerAttribute("NH");

				if(sr.getReadNegativeStrandFlag() == isNegativeStrand && NH <= max_hits){
					int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentStart() : sr.getAlignmentEnd();
					if(terminus >= start && terminus <= end){
						n++;
					}
				}

			}

			sri.close();

			return n;
		}
		if(strandedness == Strandedness.reverse){

			while(sri.hasNext()){
				SAMRecord sr = sri.next();

				int NH = sr.getIntegerAttribute("NH");

				if(sr.getReadNegativeStrandFlag() != isNegativeStrand && NH <= max_hits){
					int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentEnd() : sr.getAlignmentStart();
					if(terminus >= start && terminus <= end){
						n++;
					}
				}
			}

			sri.close();

			return n;
		}
		else{
			sri.close();
			throw new RuntimeException("Unimplemented");
		}
	}

	public static char strand(SAMRecord sr, Strandedness strandedness) {
		if(strandedness == Strandedness.unstranded){
			return '.';
		}
		else if(strandedness == Strandedness.reverse_forward){
			if(sr.getFirstOfPairFlag() || sr.getSecondOfPairFlag()){	
				return sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag() ? '+' : '-';
			}
			else{
				return '.';
			}
		}
		else if(strandedness == Strandedness.forward){
			return sr.getReadNegativeStrandFlag() ? '-' : '+';
		}
		else if(strandedness == Strandedness.reverse){
			return sr.getReadNegativeStrandFlag() ? '+' : '-';
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static boolean isNegativeStrand(SAMRecord sr, Strandedness strandedness) {

		if(strandedness == Strandedness.reverse_forward){
			if(sr.getFirstOfPairFlag() || sr.getSecondOfPairFlag()){	
				return !(sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag());
			}
			else{
				throw new RuntimeException("Read does not have mate flag");
			}
		}
		else if(strandedness == Strandedness.forward){
			return sr.getReadNegativeStrandFlag();
		}
		else if(strandedness == Strandedness.reverse){
			return !sr.getReadNegativeStrandFlag();
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static char readStrand(SAMRecord sr, Strandedness strandedness){

		if(strandedness == Strandedness.unstranded){
			return '.';
		}
		else if(strandedness == Strandedness.reverse_forward){
			if(sr.getReadPairedFlag() && (sr.getReadNegativeStrandFlag() ^ sr.getMateNegativeStrandFlag())){	
				if((sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag())){
					return '+';
				}
				else{
					return '-';
				}
			}
			else{
				return '.';
			}
		}
		else if(strandedness == Strandedness.forward){
			if(sr.getReadNegativeStrandFlag()){
				return '-';
			}
			else{
				return '+';
			}		
		}
		else if(strandedness == Strandedness.reverse){
			if(sr.getReadNegativeStrandFlag()){
				return '+';
			}
			else{
				return '-';
			}	
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static int nAlignedReads(SamReader sfr, AnnotatedRegion ar, boolean contained){
		return nAlignedReads(sfr, ar.chr, ar.start, ar.end, contained);
	}

	public static int nAlignedReads(SamReader sfr, String chr, int start, int end, boolean contained){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			sri.next();
			n++;
		}
		sri.close();
		return n;
	}

	public static int nAlignedReads(SamReader sfr, String chr, int start, int end, boolean contained, boolean isNegativeStrand){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n=0;

		while(sri.hasNext() ){
			SAMRecord sr = sri.next();
			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
				n++;
			}
		}

		sri.close();

		return n;
	}

	public static boolean moreThanNAlignedReads(SamReader sfr, String chr, int start, int end, boolean contained, int threshold){
		SAMRecordIterator sri = sfr.query(chr, start, end, contained);
		int n = 0;
		while(sri.hasNext()){
			sri.next();
			n++;
			if(n>threshold)
				break;
		}
		sri.close();
		return n > threshold;
	}

	public static List<String> referenceSequenceNames(SAMFileHeader sfh){
		List<String> referenceSequenceNames = new ArrayList<String>();

		int i=0;
		SAMSequenceRecord ssr = sfh.getSequence(i);
		while(ssr!=null){		
			referenceSequenceNames.add(ssr.getSequenceName());			
			i++;				
			ssr=sfh.getSequence(i);
		}

		return referenceSequenceNames;
	}

	public static Map<String,Integer> referenceSequenceLengths(SAMFileHeader sfh){
		Map<String,Integer> referenceSequenceLengths = new HashMap<String,Integer>();


		int i=0;
		SAMSequenceRecord ssr = sfh.getSequence(i);
		while(ssr!=null){		
			referenceSequenceLengths.put(ssr.getSequenceName(), ssr.getSequenceLength());			
			i++;				
			ssr=sfh.getSequence(i);
		}

		return referenceSequenceLengths;

	}

	/**
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param end
	 * @param contained
	 * @param strandedness
	 * @return if the library is stranded, the positive coverage is returned in index-0 of the array, negative coverage in index-1
	 * 
	 */
	public static int[][] coverage(SamReader SAMReader, String chromosome, int start, int end, boolean contained, Strandedness strandedness){
		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome, start, end, contained);
		//		}
		int length = 1+end-start;

		if(strandedness == Strandedness.unstranded){
			int[][] cov = new int[1][length];

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int alignmentPosition = sr.getAlignmentStart();
				Cigar cigar = sr.getCigar();

				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
						for(int j=0; j<cigarElement.getLength(); j++){
							if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
								cov[0][alignmentPosition-start]++;
							}
							alignmentPosition++;
						}
					}			
				}

			}

			sri.close();

			return cov;
		}
		else if(strandedness == Strandedness.reverse_forward){
			int[][] cov = new int[2][length];

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();
				if(sr.getReadPairedFlag()){
					int mate_index = sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag() ? 0 : 1;

					int alignmentPosition = sr.getAlignmentStart();
					Cigar cigar = sr.getCigar();

					for(int i=0; i<cigar.numCigarElements(); i++){
						CigarElement cigarElement = cigar.getCigarElement(i);
						if(cigarElement.getOperator().consumesReferenceBases()){
							boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
							for(int j=0; j<cigarElement.getLength(); j++){
								if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
									cov[mate_index][alignmentPosition-start]++;
								}
								alignmentPosition++;
							}
						}			
					}
				}
			}

			sri.close();

			return cov;
		}
		else if(strandedness == Strandedness.forward){
			int[][] cov = new int[2][length];

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int mate_index = sr.getReadNegativeStrandFlag() ? 1 : 0;

				int alignmentPosition = sr.getAlignmentStart();
				Cigar cigar = sr.getCigar();

				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
						for(int j=0; j<cigarElement.getLength(); j++){
							if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
								cov[mate_index][alignmentPosition-start]++;
							}
							alignmentPosition++;
						}
					}			
				}

			}

			sri.close();

			return cov;
		}
		else if(strandedness == Strandedness.reverse){
			int[][] cov = new int[2][length];

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int mate_index = sr.getReadNegativeStrandFlag() ? 0 : 1;

				int alignmentPosition = sr.getAlignmentStart();
				Cigar cigar = sr.getCigar();

				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
						for(int j=0; j<cigarElement.getLength(); j++){
							if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
								cov[mate_index][alignmentPosition-start]++;
							}
							alignmentPosition++;
						}
					}			
				}

			}

			sri.close();

			return cov;
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static int[] coverage(SamReader SAMReader, String chromosome, int start, int end, boolean contained, boolean isNegativeStrand){
		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);
		if(!sri.hasNext()){
			sri.close();
			sri = SAMReader.query(chromosome, start, end, contained);
		}
		int length = 1+end-start;
		int[] cov = new int[length];

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
			if((sr.getAttribute("XS")==null && !(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)) || (sr.getAttribute("XS")!=null && !(sr.getAttribute("XS").equals('-') ^ isNegativeStrand))){

				int alignmentPosition = sr.getAlignmentStart();
				Cigar cigar = sr.getCigar();

				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
						for(int j=0; j<cigarElement.getLength(); j++){
							if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
								cov[alignmentPosition-start]++;
							}
							alignmentPosition++;
						}
					}			
				}

				/*	for(int i=sr.getAlignmentStart(); i<=sr.getAlignmentEnd(); i++){
						if( i >= start && i <=  end){
							cov[i-start]++;
						}
					}*/
			}
		}

		sri.close();

		return cov;
	}

	/**
	 * Could round end up to the next bin that contains all the reads
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param end
	 * @param binSize
	 * @param contained
	 * @return
	 */
	public static int[] coverage(SamReader SAMReader, String chromosome, int start, int end, int binSize, boolean contained, boolean isNegativeStrand){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = 1+end-start;
		int[] cov = new int[(length+binSize-1)/binSize];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();


			if((sr.getAttribute("XS")==null && !(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)) || (sr.getAttribute("XS")!=null && !(sr.getAttribute("XS").equals('-') ^ isNegativeStrand))){

				//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){

				int alignmentPosition = sr.getAlignmentStart();
				Cigar cigar = sr.getCigar();

				for(int i=0; i<cigar.numCigarElements(); i++){
					CigarElement cigarElement = cigar.getCigarElement(i);
					if(cigarElement.getOperator().consumesReferenceBases()){
						boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
						int lastIncremented=-1;
						for(int j=0; j<cigarElement.getLength(); j++){
							if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end && lastIncremented!=((alignmentPosition-start)/binSize)){
								lastIncremented=(alignmentPosition-start)/binSize;
								cov[(alignmentPosition-start)/binSize]++;
							}
							alignmentPosition++;
						}
					}			
				}
			}
		}

		sri.close();

		return cov;
	}

	public static int[] coverage(SamReader SAMReader, String chromosome, int start, int end, int binSize, boolean contained){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = 1+end-start;
		int[] cov = new int[(length+binSize-1)/binSize];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){

			int alignmentPosition = sr.getAlignmentStart();
			Cigar cigar = sr.getCigar();

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
					int lastIncremented=-1;
					for(int j=0; j<cigarElement.getLength(); j++){
						if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end && lastIncremented!=((alignmentPosition-start)/binSize)){
							lastIncremented=(alignmentPosition-start)/binSize;
							cov[(alignmentPosition-start)/binSize]++;
						}
						alignmentPosition++;
					}
				}			
			}
		}

		sri.close();

		return cov;
	}

	public static int[] binnedMaxCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained){
		int[] cov = coverage(SAMReader, chromosome, start, start + (nBins*binSize) - 1, contained);
		int[] binnedMaxCov = new int[nBins];

		for(int i=0; i<nBins; i++){
			for(int j=0; j<binSize; j++){
				if(cov[i*binSize + j] > binnedMaxCov[i])
					binnedMaxCov[i] = cov[i*binSize + j];
			}
		}

		return binnedMaxCov;
	}

	public static int[] maxBinnedCoverage(SamReader SAMReader, String chromosome, int start, int end, int nBins, boolean contained, boolean isNegativeStrand){
		int binSize = (int) Math.ceil((end-start)/((double) nBins));
		return binnedMaxCoverage(SAMReader, chromosome, start, nBins, binSize, contained, isNegativeStrand);
	}

	public static int[] maxBinnedCoverage(SamReader SAMReader, String chromosome, int start, int end, int nBins, boolean contained){
		int binSize = (int) Math.ceil((end-start)/((double) nBins));
		return binnedMaxCoverage(SAMReader, chromosome, start, nBins, binSize, contained);
	}

	public static int[] binnedMaxCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, boolean isNegativeStrand){
		int[] cov = coverage(SAMReader, chromosome, start, start + (nBins*binSize) - 1, contained, isNegativeStrand);
		int[] binnedMaxCov = new int[nBins];

		for(int i=0; i<nBins; i++){
			for(int j=0; j<binSize; j++){
				if(cov[i*binSize + j] > binnedMaxCov[i])
					binnedMaxCov[i] = cov[i*binSize + j];
			}
		}

		return binnedMaxCov;
	}

	public static int[] binnedCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = nBins*binSize;
		int end = start+length-1;
		int[] cov = new int[nBins];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){

			int alignmentPosition = sr.getAlignmentStart();
			Cigar cigar = sr.getCigar();

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
					int lastIncremented=-1;
					for(int j=0; j<cigarElement.getLength(); j++){
						if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end && lastIncremented!=((alignmentPosition-start)/binSize)){
							lastIncremented=(alignmentPosition-start)/binSize;
							cov[(alignmentPosition-start)/binSize]++;
						}
						alignmentPosition++;
					}
				}			
			}
		}

		sri.close();

		return cov;
	}

	public static int[] binnedCoverage(Collection<SamReader> SAMReaders, String chromosome, int start, int nBins, int binSize, boolean contained){
		int[] cov = new int[nBins];

		for(SamReader sfr : SAMReaders){
			int[] tCov = binnedCoverage(sfr, chromosome, start, nBins, binSize, contained);
			for(int i=0; i<tCov.length; i++){
				cov[i] += tCov[i];
			}
		}

		return cov;
	}

	public static int[] binnedCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, boolean alignmentStart){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = nBins*binSize;
		int[] cov = new int[nBins];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
			int pos = -1;
			if(alignmentStart){
				pos = (sr.getAlignmentStart()-start)/binSize;		
			}
			else{
				pos = (sr.getAlignmentEnd()-start)/binSize;
			}

			if(pos>-1 && pos < nBins)
				cov[pos]++;


		}

		sri.close();

		return cov;
	}

	public static int[] binnedStartCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, boolean alignmentStart){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = nBins*binSize;
		int[] cov = new int[nBins];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
			int pos = -1;
			if(alignmentStart ^ sr.getReadNegativeStrandFlag()){
				pos = (sr.getAlignmentStart()-start)/binSize;		
			}
			else{
				pos = (sr.getAlignmentEnd()-start)/binSize;
			}

			if(pos>-1 && pos < nBins)
				cov[pos]++;


		}

		sri.close();

		return cov;
	}

	public static int[] binnedMaxEndCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = nBins*binSize;
		int[] start_cov = new int[nBins];
		int[] end_cov = new int[nBins];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
			int start_pos = (sr.getAlignmentStart()-start)/binSize;		
			int	end_pos = (sr.getAlignmentEnd()-start)/binSize;

			if(start_pos>-1 && start_pos < nBins)
				start_cov[start_pos]++;
			if(end_pos>-1 && end_pos < nBins)
				end_cov[end_pos]++;


		}
		sri.close();

		for(int i=0; i<nBins; i++){
			if(start_cov[i] > end_cov[i])
				end_cov[i] = start_cov[i];
		}


		return end_cov;
	}

	/**
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param nBins
	 * @param binSize
	 * @param contained
	 * @param strandedness
	 * @return if the data is stranded the positive data is in the first position of the array, negative is in the second
	 */
	public static int[] binnedMaxEndCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, Strandedness strandedness, boolean isNegativeStrand){
		if(strandedness == Strandedness.unstranded){
			return binnedMaxEndCoverage(SAMReader, chromosome, start, nBins, binSize, contained, strandedness)[0];
		}
		else{
			return binnedMaxEndCoverage(SAMReader, chromosome, start, nBins, binSize, contained, strandedness)[isNegativeStrand?1:0];
		}
	}

	public static int[] binnedMaxEndCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, boolean unique, Strandedness strandedness, boolean isNegativeStrand){
		if(strandedness == Strandedness.unstranded){
			return binnedMaxEndCoverage(SAMReader, chromosome, start, nBins, binSize, contained, unique, strandedness)[0];
		}
		else{
			return binnedMaxEndCoverage(SAMReader, chromosome, start, nBins, binSize, contained, unique, strandedness)[isNegativeStrand?1:0];
		}
	}
	/**
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param nBins
	 * @param binSize
	 * @param contained
	 * @param strandedness
	 * @return if the data is stranded the positive data is in the first position of the array, negative is in the second
	 */
	public static int[][] binnedMaxEndCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, Strandedness strandedness){
		int length = nBins*binSize;

		//		if(strandedness == Strandedness.unstranded){
		//
		//
		//		}
		//		else 
		if(strandedness == Strandedness.reverse_forward){
			int[][] start_cov = new int[2][nBins];
			int[][] end_cov = new int[2][nBins];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				if(sr.getReadPairedFlag()){
					int mate_index = sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag() ? 0 : 1;

					int start_pos = (sr.getAlignmentStart()-start)/binSize;		
					int	end_pos = (sr.getAlignmentEnd()-start)/binSize;

					if(start_pos>-1 && start_pos < nBins)
						start_cov[mate_index][start_pos]++;
					if(end_pos>-1 && end_pos < nBins)
						end_cov[mate_index][end_pos]++;
				}
			}
			sri.close();

			for(int i=0; i<nBins; i++){
				for(int j=0; j<2; j++){
					if(start_cov[j][i] > end_cov[j][i])
						end_cov[j][i] = start_cov[j][i];
				}
			}

			return end_cov;
		}
		if(strandedness == Strandedness.unstranded){
			int[] start_cov = new int[nBins];
			int[] end_cov = new int[nBins];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int start_pos = (sr.getAlignmentStart()-start)/binSize;		
				int	end_pos = (sr.getAlignmentEnd()-start)/binSize;

				if(start_pos>-1 && start_pos < nBins)
					start_cov[start_pos]++;
				if(end_pos>-1 && end_pos < nBins)
					end_cov[end_pos]++;

			}
			sri.close();

			int[][] max_cov = new int[1][nBins];
			for(int i=0; i<nBins; i++){
				if(start_cov[i] > end_cov[i])
					max_cov[0][i] = start_cov[i];
				else
					max_cov[0][i] = end_cov[i];

			}

			return max_cov;
		}
		//		else if(strandedness == Strandedness.forward){
		//
		//
		//
		//		}
		//		else if(strandedness == Strandedness.reverse){
		//
		//
		//
		//		}

		throw new RuntimeException("Unimplemented");
	}

	public static int[][] binnedMaxEndCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, boolean unique, Strandedness strandedness){
		int length = nBins*binSize;

		//		if(strandedness == Strandedness.unstranded){
		//
		//
		//		}
		//		else 
		if(strandedness == Strandedness.reverse_forward){
			int[][] start_cov = new int[2][nBins];
			int[][] end_cov = new int[2][nBins];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();
				Integer NH = sr.getIntegerAttribute("NH");
				if(sr.getReadPairedFlag()){
					if(!unique || NH==null || NH<3){
						int mate_index = sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag() ? 0 : 1;

						int start_pos = (sr.getAlignmentStart()-start)/binSize;		
						int	end_pos = (sr.getAlignmentEnd()-start)/binSize;

						if(start_pos>-1 && start_pos < nBins)
							start_cov[mate_index][start_pos]++;
						if(end_pos>-1 && end_pos < nBins)
							end_cov[mate_index][end_pos]++;
					}
				}
			}
			sri.close();

			for(int i=0; i<nBins; i++){
				for(int j=0; j<2; j++){
					if(start_cov[j][i] > end_cov[j][i])
						end_cov[j][i] = start_cov[j][i];
				}
			}

			return end_cov;
		}
		if(strandedness == Strandedness.unstranded){
			int[] start_cov = new int[nBins];
			int[] end_cov = new int[nBins];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();
				Integer NH = sr.getIntegerAttribute("NH");
				if(!unique || NH==null || NH==1){
					int start_pos = (sr.getAlignmentStart()-start)/binSize;		
					int	end_pos = (sr.getAlignmentEnd()-start)/binSize;

					if(start_pos>-1 && start_pos < nBins)
						start_cov[start_pos]++;
					if(end_pos>-1 && end_pos < nBins)
						end_cov[end_pos]++;
				}
			}
			sri.close();

			int[][] max_cov = new int[1][nBins];
			for(int i=0; i<nBins; i++){
				if(start_cov[i] > end_cov[i])
					max_cov[0][i] = start_cov[i];
				else
					max_cov[0][i] = end_cov[i];

			}

			return max_cov;
		}
		//		else if(strandedness == Strandedness.forward){
		//
		//
		//
		//		}
		//		else if(strandedness == Strandedness.reverse){
		//
		//
		//
		//		}

		throw new RuntimeException("Unimplemented");
	}

	/**
	 * 
	 * @param sr
	 * @param strandedness
	 * @param end5p - the end in transcript coordinate space
	 * @return
	 */
	public static int terminalAlignedPosition(SAMRecord sr, Strandedness strandedness, boolean end5p){
		if(strandedness == Strandedness.forward){
			int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentStart() : sr.getAlignmentEnd();
			return terminus;
		}
		else if(strandedness == Strandedness.reverse){
			int terminus = sr.getReadNegativeStrandFlag() ^ end5p ? sr.getAlignmentEnd() : sr.getAlignmentStart();
			return terminus;
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static int[][] endCoverage(SamReader SAMReader, String chromosome, int start, int end, Strandedness strandedness, boolean contained, boolean alignmentStart){

		if(strandedness==Strandedness.reverse){
			int length = end-start+1;
			int[][] cov = new int[2][length];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int mate_index = sr.getReadNegativeStrandFlag() ? 0 : 1;


				int pos = -1;
				if(alignmentStart){
					pos = (sr.getAlignmentStart()-start);		
				}
				else{
					pos = (sr.getAlignmentEnd()-start);
				}

				if(pos>-1 && pos < length)
					cov[mate_index][pos]++;

			}

			sri.close();

			return cov;
		}
		if(strandedness==Strandedness.forward){
			int length = end-start+1;
			int[][] cov = new int[2][length];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int mate_index = sr.getReadNegativeStrandFlag() ? 1 : 0;


				int pos = -1;
				if(alignmentStart){
					pos = (sr.getAlignmentStart()-start);		
				}
				else{
					pos = (sr.getAlignmentEnd()-start);
				}

				if(pos>-1 && pos < length)
					cov[mate_index][pos]++;

			}

			sri.close();

			return cov;
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static int[][] terminalCoverage(SamReader SAMReader, String chromosome, int start, int end, Strandedness strandedness, boolean contained, boolean end5p){

		if(strandedness==Strandedness.reverse){
			int length = end-start+1;
			int[][] cov = new int[2][length];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int mate_index = sr.getReadNegativeStrandFlag() ? 0 : 1;

				int pos = terminalAlignedPosition(sr, strandedness, end5p);

				if(pos>=start && pos <= end)
					cov[mate_index][pos-start]++;

			}

			sri.close();

			return cov;
		}
		if(strandedness==Strandedness.forward){
			int length = end-start+1;
			int[][] cov = new int[2][length];

			SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);

			while(sri.hasNext() ){
				SAMRecord sr = sri.next();

				int mate_index = sr.getReadNegativeStrandFlag() ? 1 : 0;

				int pos = terminalAlignedPosition(sr, strandedness, end5p);

				if(pos>=start && pos <= end)
					cov[mate_index][pos-start]++;

			}

			sri.close();

			return cov;
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

	public static int[] binnedEndCoverage(SamReader SAMReader, String chromosome, int start, int nBins, int binSize, boolean contained, boolean alignmentStart){
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = nBins*binSize;
		int[] cov = new int[nBins];

		SAMRecordIterator sri = SAMReader.query(chromosome, start, start+length-1, contained);

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
			int pos = -1;
			if(alignmentStart){
				pos = (sr.getAlignmentStart()-start)/binSize;		
			}
			else{
				pos = (sr.getAlignmentEnd()-start)/binSize;
			}

			if(pos>-1 && pos < nBins)
				cov[pos]++;


		}

		sri.close();

		return cov;
	}

	/**
	 * 
	 * @param sr
	 * @param position
	 * @return the next reference position consumed by sr after position
	 * or -1 if the read doesn't consume any reference base after position
	 */
	public static int nextConsumedReferencePosition(SAMRecord sr, int position){

		int alignmentPosition = sr.getAlignmentStart();
		Cigar cigar = sr.getCigar();

		for(int i=0; i<cigar.numCigarElements(); i++){
			CigarElement cigarElement = cigar.getCigarElement(i);
			if(cigarElement.getOperator().consumesReferenceBases()){
				boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();

				for(int j=0; j<cigarElement.getLength(); j++){
					if(consumesReadBases && alignmentPosition>position)
						return alignmentPosition;
					alignmentPosition++;
				}
			}			
		}

		return -1;
	}

	public static int[] coverage(SamReader SAMReader, String chromosome, int start, int end, boolean contained){
		return coverage(SAMReader, chromosome, start, end, contained, false, false, true);
	}

	/**
	 * 
	 * TODO doesn't handle the isStranded option, should ignore isNegative strand if false
	 * 
	 * @param SAMReader
	 * @param chromosome
	 * @param start
	 * @param end
	 * @param contained
	 * @param isStranded
	 * @param isNegativeStrand
	 * @param removeDuplicates
	 * @return
	 */

	public static int[] coverage(SamReader SAMReader, String chromosome, int start, int end, boolean contained, boolean isStranded, boolean isNegativeStrand, boolean removeDuplicates){
		SAMRecordIterator sri = SAMReader.query(chromosome, start, end, contained);
		//		if(!sri.hasNext()){
		//			sri.close();
		//			sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		//		}
		int length = 1+end-start;
		int[] cov = new int[length];

		//		System.out.println(chromosome+":"+(center-nBasePairs)+"-"+(center+nBasePairs));
		while(sri.hasNext() ){
			SAMRecord sr = sri.next();

			if(removeDuplicates && sr.getDuplicateReadFlag()){
				//				System.out.printf("%s:%d-%d\n", sr.getReferenceName(),sr.getAlignmentStart(),sr.getAlignmentEnd());
				continue;
			}

			//			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){

			int alignmentPosition = sr.getAlignmentStart();
			Cigar cigar = sr.getCigar();

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
					for(int j=0; j<cigarElement.getLength(); j++){
						if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
							cov[alignmentPosition-start]++;
						}
						alignmentPosition++;
					}
				}			
			}

			/*	for(int i=sr.getAlignmentStart(); i<=sr.getAlignmentEnd(); i++){
						if( i >= start && i <=  end){
							cov[i-start]++;
						}
					}*/
		}

		sri.close();

		return cov;
	}

	/**
	 * 
	 * @param coverage
	 * @param depth
	 * @return fraction of the coverage that meets at least the specified depth
	 * @throws IOException
	 */
	public static double fractionCovered(int[] coverage, int depth) throws IOException{
		int n=0;

		int length=coverage.length;

		for(int i=0;i<coverage.length;i++){
			if(coverage[i]>=depth)
				n++;
		}

		return n*1./length;
	}

	public static int getNReads(SamReader SAMReader, String chromosome, int start, int end, boolean contained){
		SAMRecordIterator sri = SAMReader.query(chromosome.substring(3), start, end, contained);

		int nReads=0;

		while(sri.hasNext() ){
			sri.next();
			nReads++;
		}
		sri.close();

		return nReads;
	}

	public static int getNReads(SamReader SAMReader, String chromosome, int start, int end, boolean contained, boolean isNegativeStrand){
		SAMRecordIterator sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		int nReads=0;

		while(sri.hasNext() ){
			SAMRecord sr = sri.next();
			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
				nReads++;
			}
		}
		sri.close();
		return nReads;
	}

	public static Iterator<SAMRecord> query(SamReader SAMReader, String chromosome, int start, int end, boolean contained, boolean isNegativeStrand){
		SAMRecordIterator sri = SAMReader.query(chromosome.substring(3), start, end, contained);
		int nReads=0;

		while(sri.hasNext() ){
			SAMRecord sr = sri.next();
			if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
				nReads++;
			}
		}
		sri.close();
		return null; 
	}

	public static class StrandedSAMRecordIterator implements Iterator<SAMRecord>{
		boolean isNegativeStrand;
		SAMRecordIterator sri;
		SAMRecord next;

		public StrandedSAMRecordIterator(SamReader SAMReader, String chromosome, int start, int end, boolean contained, boolean isNegativeStrand) {
			sri = SAMReader.query(chromosome, start, end, contained);
			next = prepareNext();
		}

		public boolean hasNext() {
			return next != null;
		}

		public SAMRecord next() {
			SAMRecord toReturn = next;
			next = prepareNext();
			return toReturn;
		}

		public void remove() {

		}

		public SAMRecord prepareNext(){
			SAMRecord next = null;
			while(sri.hasNext()){
				SAMRecord sr = sri.next();
				if(!(sr.getReadNegativeStrandFlag() ^ isNegativeStrand)){
					next = sr;
					break;
				}
			}
			return next;
		}

		public void close(){
			sri.close();
		}
	}



}