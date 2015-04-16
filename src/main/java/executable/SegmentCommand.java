package executable;

import java.io.File;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.converters.BooleanConverter;
import com.beust.jcommander.converters.DoubleConverter;
import com.beust.jcommander.converters.IntegerConverter;

@XmlRootElement(name="SegmentCommand")
@XmlAccessorType(XmlAccessType.PROPERTY)
public class SegmentCommand{	
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
	
	@XmlElement
	public void setBam(File bam) {
		this.bam = bam;
	}

	@XmlElement
	public void setStrandedness(String strandedness) {
		this.strandedness = strandedness;
	}

	@XmlElement
	public void setDir(String dir) {
		this.dir = dir;
	}

	@XmlElement
	public void setBase(String base) {
		this.base = base;
	}

	@XmlElement
	public void setW(int w) {
		this.w = w;
	}

	@XmlElement
	public void setSegment_r(int segment_r) {
		this.segment_r = segment_r;
	}

	@XmlElement
	public void setSegment_p(double segment_p) {
		this.segment_p = segment_p;
	}

	@XmlElement
	public void setNb_r(int nb_r) {
		this.nb_r = nb_r;
	}

	@XmlElement
	public void setMerge_radius(int merge_radius) {
		this.merge_radius = merge_radius;
	}

	@XmlElement
	public void setInternal(boolean internal) {
		this.internal = internal;
	}

	@XmlElement
	public void setCoverage(boolean coverage) {
		this.coverage = coverage;
	}

	@XmlElement
	public void setMin_fold(double min_fold) {
		this.min_fold = min_fold;
	}

	@XmlElement
	public void setMin_terminal(int min_terminal) {
		this.min_terminal = min_terminal;
	}

	@XmlElement
	public void setJnct_alpha(double jnct_alpha) {
		this.jnct_alpha = jnct_alpha;
	}

	@XmlElement
	public void setFilled_gap_segments(String filled_gap_segments) {
		this.filled_gap_segments = filled_gap_segments;
	}

	public File getBam() {
		return bam;
	}

	public String getStrandedness() {
		return strandedness;
	}

	public String getDir() {
		return dir;
	}

	public String getBase() {
		return base;
	}

	public int getW() {
		return w;
	}

	public int getSegment_r() {
		return segment_r;
	}

	public double getSegment_p() {
		return segment_p;
	}

	public int getNb_r() {
		return nb_r;
	}

	public int getMerge_radius() {
		return merge_radius;
	}

	public boolean isInternal() {
		return internal;
	}

	public boolean isCoverage() {
		return coverage;
	}

	public double getMin_fold() {
		return min_fold;
	}

	public int getMin_terminal() {
		return min_terminal;
	}

	public double getJnct_alpha() {
		return jnct_alpha;
	}

	public String getFilled_gap_segments() {
		return filled_gap_segments;
	}
}
