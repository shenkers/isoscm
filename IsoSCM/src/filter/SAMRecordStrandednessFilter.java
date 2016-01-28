package filter;

import htsjdk.samtools.SAMRecord;
import tools.Strandedness;

public class SAMRecordStrandednessFilter extends ComposableFilter<SAMRecord>{

	Strandedness strandedness;
	public boolean isNegativeStrand;

	public SAMRecordStrandednessFilter(Strandedness strandedness, boolean isNegativeStrand) {
		super();
		this.strandedness = strandedness;
		this.isNegativeStrand = isNegativeStrand;
	}

	public SAMRecordStrandednessFilter(Strandedness strandedness) {
		super();
		this.strandedness = strandedness;
	}

	public boolean satisfies(SAMRecord sr) {
		if(strandedness == Strandedness.unstranded){
			return true;
		}
		else if(strandedness == Strandedness.reverse_forward){
			return (sr.getReadNegativeStrandFlag() ^ sr.getSecondOfPairFlag()) != isNegativeStrand;
		}
		else if(strandedness == Strandedness.forward){
			return sr.getReadNegativeStrandFlag() == isNegativeStrand;
		}
		else if(strandedness == Strandedness.reverse){
			return sr.getReadNegativeStrandFlag() != isNegativeStrand;
		}
		else{
			throw new RuntimeException("Unimplemented");
		}
	}

}
