package processing;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import tools.AnnotatedRegion;
import tools.BAMTools;
import tools.BEDTools.BEDWriter;
import tools.StrandedGenomicIntervalSet;
import tools.Strandedness;
import util.Util;
import util.Util.ExtremeTracker;

public class SlidingWindow {

	/**
	 * Takes a bam file, and outputs all genomic segments that meet a minimum coverage criteria.
	 * The criteria is defined by a window size, theshold number of reads, and number of positions
	 * in the window that meet this criteria.
	 * 
	 * @param sfr
	 * @param w 
	 * @param thresholds 
	 * @param positionsMeetingThreshold 
	 * @param baseName
	 * @param isNegativeStrand 
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static void identifyExpressed(SAMFileReader sfr, Strandedness strandedness, int positionsMeetingThreshold, int threshold, int w, BEDWriter bedWriter) throws FileNotFoundException, IOException{

		//		String chr = "chr3";
		//		int start = 144085495;
		//		int end = 144098053;


		// window size
		//		int w=100;
		// number of positions meeting threshold
		//		int positionsMeetingThreshold=80;
		// can't read all the reads into memory at once, this is the chunksize read
		int chunkSize = 2000;

		// since reading the bam file is the bottleneck, you can output segments for different thresholds simultaneously 
		//		List<Integer> thresholds = Util.list(1,5,10,20,40,60,100);

		// get the length of each chromosome from the bamfile
		Map<String,Integer> refLengths = BAMTools.referenceSequenceLengths(sfr.getFileHeader());

		// holds the positions that satisfy the minimum expression criteria in the previous w positions
		boolean[][] areSatisfied = null;
		if(strandedness==Strandedness.unstranded){
			areSatisfied = new boolean[1][w];
		}
		else if(strandedness==Strandedness.reverse_forward){
			areSatisfied = new boolean[2][w];
		}
		// holds the segments that meet the criteria
		StrandedGenomicIntervalSet satisfiedIntervals = new StrandedGenomicIntervalSet();

		for(String chr : refLengths.keySet()){
			//			for(String chr : Util.list("3")){
			//					System.out.println(chr);
			identifyChrExpressed(chr, sfr, strandedness, chunkSize, w, positionsMeetingThreshold, threshold, satisfiedIntervals, bedWriter);			
		}

	}

	public static void identifyChrExpressed(String chr, SAMFileReader sfr, Strandedness strandedness, int chunkSize, int w, int positionsMeetingThreshold, int threshold, StrandedGenomicIntervalSet satisfiedIntervals, BEDWriter bedWriter){

		Map<String,Integer> refLengths = BAMTools.referenceSequenceLengths(sfr.getFileHeader());
		int refLength = refLengths.get(chr);

		// get the next read, this lets us skip all chunks that have no reads
		Integer nextStart = null;
		{
			SAMRecordIterator sri = sfr.queryOverlapping(chr, 1, refLength);

			if(sri.hasNext()){
				SAMRecord next = sri.next();
				nextStart = next.getAlignmentStart();
			}

			sri.close();
		}

		//		System.out.println(nextStart);
		if(nextStart!=null){
			for(int chunkStart=nextStart; chunkStart<refLength; chunkStart+=chunkSize-w){

				// if there are any segments that are more than w away from the chunk start they can be written


				/* TODO need to fix remove (SEE NOTEBOOK) for this to work*/
				List<AnnotatedRegion> written = new LinkedList<AnnotatedRegion>();
				
				for(AnnotatedRegion r : satisfiedIntervals){
					if(r.end < chunkStart-w){
						written.add(r);
						//						System.out.println(written);
						bedWriter.write(r.chr, r.start + (w-positionsMeetingThreshold), r.end - (w-positionsMeetingThreshold), r.strand);
					}
				}

				for(AnnotatedRegion r : written){
					satisfiedIntervals.remove(chr, r.start, r.end, r.strand);
				}
				/* */

				int chunkEnd = Math.max(Math.min(chunkStart+chunkSize-1, refLength),chunkStart+w);

				Iterator<AnnotatedRegion> it = satisfiedIntervals.iterator();
				while(it.hasNext()){
					it.next();
				}
				
//								bedWriter.write(chr, chunkStart, chunkEnd, '.');
				// if there are no reads in this chunk, skip it
				//				if(nextStart != null && nextStart > chunkEnd){
				//					//					System.out.println("continuing "+chunkStart);
				//					continue;
				//				}

				boolean alignedReads = false;

				// calculate the coverage in this chunk
				int[][] cov = BAMTools.coverage(sfr, chr, chunkStart, chunkEnd, false, strandedness);

				if(strandedness==Strandedness.unstranded){
					// keeps a cumulative sum of positions that meet the criteria in the previous w positions
					int nSatisfied = 0;
					// TODO clear array rather than create new
					boolean[] areSatisfied = new boolean[w];

					// check which positions satisfy the threshold in the first w positions
					for(int i=0; i<w-1; i++){
						boolean isPosSatisfied = cov[0][i]>=threshold;
						
						alignedReads |= (cov[0][i]>0);

						areSatisfied[i%w] = isPosSatisfied;
					
						if(isPosSatisfied)
							nSatisfied++;
					}


					// for position w+1 to the end of the chunk, add all intervals that satisfy the threshold
					for (int i = w-1; i < chunkEnd - chunkStart + 1 ; i++) {
						boolean isSatisfied = cov[0][i]>=threshold;
					
						alignedReads |= (cov[0][i]>0);

						if(areSatisfied[i%w]^isSatisfied){
							if(isSatisfied)
								nSatisfied++;
							else
								nSatisfied--;
						}
						
						areSatisfied[i%w]=isSatisfied;
						
						if(nSatisfied >= positionsMeetingThreshold){
							satisfiedIntervals.add(chr, chunkStart+i-w+1, chunkStart + i, '.');
						}
					}
				}
				else if(strandedness==Strandedness.reverse_forward || strandedness==Strandedness.forward || strandedness == Strandedness.reverse){
					// keeps a cumulative sum of positions that meet the criteria in the previous w positions
					int nPosSatisfied = 0;
					int nNegSatisfied = 0;
					// TODO clear array rather than create new
					boolean[][] areSatisfied = new boolean[2][w];

					// check which positions satisfy the threshold in the first w positions
					for(int i=0; i<w-1; i++){
						boolean isPosSatisfied = cov[0][i]>=threshold;
						boolean isNegSatisfied = cov[1][i]>=threshold;

						alignedReads |= (cov[0][i]>0||cov[1][i]>0);

						areSatisfied[0][i%w] = isPosSatisfied;
						areSatisfied[1][i%w] = isNegSatisfied;

						if(isPosSatisfied)
							nPosSatisfied++;
						if(isNegSatisfied)
							nNegSatisfied++;
					}


					// for position w+1 to the end of the chunk, add all intervals that satisfy the threshold
					for (int i = w-1; i < chunkEnd - chunkStart + 1 ; i++) {
						boolean isPosSatisfied = cov[0][i]>=threshold;
						boolean isNegSatisfied = cov[1][i]>=threshold;

						alignedReads |= (cov[0][i]>0||cov[1][i]>0);

						//			System.out.println(nS);
						if(areSatisfied[0][i%w]^isPosSatisfied){
							if(isPosSatisfied)
								nPosSatisfied++;
							else
								nPosSatisfied--;
						}
						if(areSatisfied[1][i%w]^isNegSatisfied){
							if(isNegSatisfied)
								nNegSatisfied++;
							else
								nNegSatisfied--;
						}

						areSatisfied[0][i%w]=isPosSatisfied;
						areSatisfied[1][i%w]=isNegSatisfied;

						//						if(isNegSatisfied)
						//							bedWriter.write(chr, chunkStart+i-w+1, chunkStart + i, '-');

						if(nPosSatisfied >= positionsMeetingThreshold){
							//							System.out.printf("adding %d\n", chunkStart+i-w+1);
							//							System.out.println(satisfiedIntervals.size());
							satisfiedIntervals.add(chr, chunkStart+i-w+1, chunkStart + i, '+');
							//							bedWriter.write(chr, chunkStart+i-w+1, chunkStart + i, '+');
						}
						if(nNegSatisfied >= positionsMeetingThreshold){
							//							System.out.printf("adding %d\n", chunkStart+i-w+1);
							//							System.out.println(satisfiedIntervals.size());
							satisfiedIntervals.add(chr, chunkStart+i-w+1, chunkStart + i, '-');
							//							bedWriter.write(chr, chunkStart+i-w+1, chunkStart + i, '-');

						}
					}
				}
				else{
					throw new RuntimeException("strandedness "+strandedness+" not handled.");
				}

				// check where the next read is so we can skip empty chunks
				if(!alignedReads && chunkEnd < refLength){
					ExtremeTracker<Integer> nextConsumedBase = new ExtremeTracker<Integer>(new Util.ComparableComparator<Integer>());

					SAMRecordIterator sri = sfr.queryContained(chr, chunkEnd, refLength);

					if(sri.hasNext()){
						SAMRecord next = sri.next();
						nextStart = next.getAlignmentStart();
						nextConsumedBase.put(nextStart);
					}

					sri.close();

					sri = sfr.queryOverlapping(chr, chunkEnd, chunkEnd);
					while(sri.hasNext()){
						SAMRecord next = sri.next();
						int nextConsumed = BAMTools.nextConsumedReferencePosition(next, chunkEnd);
						if(nextConsumed!=-1){
							nextConsumedBase.put(nextConsumed);
						}
					}
					sri.close();

					if(nextConsumedBase.hasExtrema){
//						System.out.printf("nextstart \t %d\n", nextConsumedBase.getMin());
						chunkStart=nextConsumedBase.getMin()-chunkSize;
					}
					else{
//						System.out.printf("breaking %d \n", chunkStart);
						break;
					}

				}
			}
		}

		//		System.out.println(gis.size());
		//		Util.save(gis, fileName+".obj");		
		//		System.out.println(satisfiedIntervals.size());
		//		System.out.println(satisfiedIntervals.iterator().hasNext());
		//		System.out.println(satisfiedIntervals.iterator().next());
		// write all unwritten segments
		List<AnnotatedRegion> written = new LinkedList<AnnotatedRegion>();
		for(AnnotatedRegion r : satisfiedIntervals){
			written.add(r);
			bedWriter.write(r.chr, r.start + (w-positionsMeetingThreshold), r.end - (w-positionsMeetingThreshold), r.strand);
		}

		for(AnnotatedRegion r : written){
			satisfiedIntervals.remove(chr, r.start, r.end, r.strand);
		}

	}

	public static void main(String[] args) throws FileNotFoundException, IOException {
		{
		SAMFileReader sfr = new SAMFileReader(new File("/home/sol/lailab/sol/sangercenter/de_dup/hippocampus.bam"));
		BEDWriter bw = new BEDWriter("unstranded.bed");
		StrandedGenomicIntervalSet satisfiedIntervals = new StrandedGenomicIntervalSet();
		int threshold=1;
		int positionsMeetingThreshold=1;
		int w=1;
		int chunkSize=10000;
		Strandedness strandedness=Strandedness.unstranded;
		identifyChrExpressed("chr6", sfr, strandedness, chunkSize, w, positionsMeetingThreshold, threshold, satisfiedIntervals, bw);
		}
		
		if(false){
		SAMFileReader sfr = new SAMFileReader(new File("/mnt/LaiLab/sol/GSE41637/mapped/indexed/SRR594398.bam"));
		BEDWriter bw = new BEDWriter("/mnt/LaiLab/sol/GSE41637/isoscm/test.bed");
		identifyExpressed(sfr, Strandedness.reverse_forward, 1, 1, 1, bw);
		bw.close();
	}
	}

}
