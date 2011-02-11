package tools.shoreMap;

public class Mutation implements Comparable<Mutation>{

	private String chr;
	private int pos;
	private char newNuc;
	
	public Mutation(String chr, int pos, char newNuc) {
		super();
		this.chr = chr;
		this.pos = pos;
		this.newNuc = newNuc;
	}
	
	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public char getNewNuc() {
		return newNuc;
	}

	public int compareTo(Mutation o) {
		if(this.chr.equals(o.chr)){
			return this.chr.compareTo(o.chr);
		}
		return this.pos-o.pos;
	}
	
	
}
