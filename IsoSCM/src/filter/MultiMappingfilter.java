package filter;

import htsjdk.samtools.SAMRecord;

/**
 * 
 * Removes all reads that map fewer than a specified
 * number of times
 *
 */
public class MultiMappingfilter extends ComposableFilter<SAMRecord> {

	int n;
	
	public MultiMappingfilter(int n) {
		super();
		this.n = n;
	}

	public boolean satisfies(SAMRecord sr) {
		return sr.getAttribute("NH") != null && sr.getIntegerAttribute("NH") <= n;
	}
}
