package filter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ConsumingReadFilter extends ComposableFilter<SAMRecord> {

	public String chr;
	public int start, end;

	public ConsumingReadFilter(String chr, int start, int end) {
		super();
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	public boolean satisfies(SAMRecord sr) {
		if(sr.getReferenceName().equals(chr)){
			int alignmentPosition = sr.getAlignmentStart();
			Cigar cigar = sr.getCigar();

			for(int i=0; i<cigar.numCigarElements(); i++){
				CigarElement cigarElement = cigar.getCigarElement(i);
				if(cigarElement.getOperator().consumesReferenceBases()){
					boolean consumesReadBases = cigarElement.getOperator().consumesReadBases();
					for(int j=0; j<cigarElement.getLength(); j++){
						if(consumesReadBases && alignmentPosition >= start && alignmentPosition <= end){
							return true;
						}
						alignmentPosition++;
					}
				}			
			}			
		}
		return false;
	}
}
