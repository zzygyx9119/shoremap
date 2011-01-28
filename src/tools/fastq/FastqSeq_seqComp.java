package tools.fastq;

import java.util.Comparator;

public class FastqSeq_seqComp implements Comparator<FastqSeq> {

	public int compare(FastqSeq arg0, FastqSeq arg1) {
		// TODO Auto-generated method stub
		return arg0.getSeq().compareToIgnoreCase(arg1.getSeq());
	}

}
