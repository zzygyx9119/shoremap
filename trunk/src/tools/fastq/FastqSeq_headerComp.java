package tools.fastq;

import java.util.Comparator;

public class FastqSeq_headerComp implements Comparator<FastqSeq> {

	public int compare(FastqSeq o1, FastqSeq o2) {
		// TODO Auto-generated method stub
		return o1.getHeader().compareTo(o2.getHeader());
	}

}
