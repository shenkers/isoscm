package tools;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import util.Util;

public class BEDTools {

	/*
	 *     chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 

	The 9 additional optional BED fields are:

    name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
    shade 	  	  	  	  	  	  	  	  	 
    score in range   	≤ 166 	167-277 	278-388 	389-499 	500-611 	612-722 	723-833 	834-944 	≥ 945
    strand - Defines the strand - either '+' or '-'.
    thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
    thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
    itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
    blockCount - The number of blocks (exons) in the BED line.
    blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
    blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
	 */
	public static class BEDWriter{
		PrintStream out;
		
		public BEDWriter(){
			out = System.out;
		}
		
		public BEDWriter(String fileName) throws FileNotFoundException {
			out = new PrintStream(new File(fileName));	
		}
		
		public BEDWriter(PrintStream out){
			this.out = out;
		}
		
		public BEDWriter(String fileName, String trackName, String description, boolean useColors) throws FileNotFoundException {
			out = new PrintStream(new File(fileName));	
			out.printf("track name=\"%s\" description=\"%s\" itemRgb=\"%s\" visibility=2\n",trackName,description,useColors?"On":"Off");
		}
		
		public void write(String chr, int start, int end, char strand){
			out.printf("%s\t%d\t%d\t\t0\t%s\n", chr,start-1,end,strand);
		}
		
		public void write(String name, String chr, int start, int end, char strand){
			out.printf("%s\t%d\t%d\t%s\t0\t%s\n", chr,start-1,end,name,strand);
		}
		
		public void write(String name, String chr, int start, int end, List<Integer> exon_starts, List<Integer> exon_lengths, char strand){
			exon_lengths = Util.sort(exon_lengths, exon_starts, true);
			Collections.sort(exon_starts);
			List<Integer> relative_starts = new ArrayList<Integer>(exon_starts.size());
			for(Integer exon_start : exon_starts){
				relative_starts.add(exon_start - start);
			}
			out.printf("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0,0,0\t%d\t%s\t%s\n", chr,start-1,end,name,strand,start-1,end,exon_starts.size(),StringUtils.join(exon_lengths,","),StringUtils.join(relative_starts,","));
		}
		
		public void write(String name, String chr, int start, int end, int thickStart, int thickEnd, List<Integer> exon_starts, List<Integer> exon_lengths, char strand){
			exon_lengths = Util.sort(exon_lengths, exon_starts, true);
			Collections.sort(exon_starts);
			List<Integer> relative_starts = new ArrayList<Integer>(exon_starts.size());
			for(Integer exon_start : exon_starts){
				relative_starts.add(exon_start - start);
			}
			out.printf("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0,0,0\t%d\t%s\t%s\n", chr,start-1,end,name,strand,thickStart-1,thickEnd,exon_starts.size(),StringUtils.join(exon_lengths,","),StringUtils.join(relative_starts,","));
		}
		
		public void write(String name, String chr, int start, int end, double score, char strand){
			out.printf("%s\t%d\t%d\t%s\t%f\t%s\n", chr,start-1,end,name, score, strand);
		}
		
		public void write(String name, String chr, int start, int end, char strand, Color color){
			out.printf("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t%s\n", chr,start-1,end,name, strand, start-1, end,StringUtils.join(Util.list(color.getRed(),color.getGreen(),color.getBlue()),","));
		}
		
		public void close(){
			out.close();
		}
	}
}
