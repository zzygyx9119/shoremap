package tools.fastq;


public class FastqCorrection implements Comparable<FastqCorrection>{

	private int pos,qual;
	private char to,from;
	
	public FastqCorrection(int pos, int qual, char to, char from) {
		super();
		this.pos = pos;
		this.qual = qual;
		this.to = to;
		this.from = from;
	}

	public int getPos() {
		return pos;
	}

	public int getQual() {
		return qual;
	}

	public char getTo() {
		return to;
	}

	public char getFrom() {
		return from;
	}
	
	public String toString(){
		return "Pos: "+getPos()+": "+getFrom()+" to "+getTo()+", qual: "+getQual();
	}

	public int compareTo(FastqCorrection arg0) {
		return this.getPos()-arg0.getPos();
	}
	
	
}
