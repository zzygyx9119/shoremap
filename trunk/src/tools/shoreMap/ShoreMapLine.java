package tools.shoreMap;

import java.util.ArrayList;
import java.util.Collections;

public class ShoreMapLine {

	private String chr;
	private int pos;
	private String alignment;
	private String readID;
	private String strand;
	private int mismatches;
	private int hits;
	private int read_length;
	private int offset;
	private String pe_flag;
	private String quality;
	
	
	public ShoreMapLine(String chr, int pos, String alignment, String readID,
			String strand, int mismatches, int hits, int read_length,
			int offset, String pe_flag, String quality) {
		super();
		this.chr = chr;
		this.pos = pos;
		this.alignment = alignment;
		this.readID = readID;
		this.strand = strand;
		this.mismatches = mismatches;
		this.hits = hits;
		this.read_length = read_length;
		this.offset = offset;
		this.pe_flag = pe_flag;
		this.quality = quality;
	}
	
	public ShoreMapLine(String[] l){
		this(l[0], Integer.parseInt(l[1]), l[2], l[3], l[4], Integer.parseInt(l[5]), Integer.parseInt(l[6]), Integer.parseInt(l[7]), Integer.parseInt(l[8]), l[9], l[10]);
	}
	
	public ShoreMapLine(String line){
		this(line.split("\t"));
	}
	
	public boolean onChr(String chr){
		return this.chr.equals(chr);
	}
	
	public boolean containsPos(int pos){
		return (pos>=this.pos&&pos<(this.pos+this.read_length));
	}
	
	public String getSeq(){
		String seq="";
		for(int i=0;i<alignment.length();i++){
			switch (alignment.charAt(i)) {
			case 'A':
			case 'T':
			case 'G':
			case 'C':
			case 'N':
				seq+=""+alignment.charAt(i);
				break;
			case '[':
				seq+=""+alignment.charAt(i+2);
				i+=3;
				break;
			default:
				System.err.println("Unknown alignment character: ("+alignment.charAt(i)+"), in read: \n"+this.toString());
				break;
			}
		}
		return seq;
	}
	
	public String getMutatedSeq(ArrayList<Mutation> mutations){
		String orig=this.getSeq();
		ArrayList<Mutation> overlapping= new ArrayList<Mutation>();
		for (Mutation mutation : mutations) {
			if(onChr(mutation.getChr())&&containsPos(mutation.getPos())){
				overlapping.add(mutation);
			}
		}
		if(overlapping.size()>0){
			Collections.sort(overlapping);
			String seq="";
			int lastPos=0;
			for (Mutation mutation : overlapping) {
				int curPos=chrPosToReadPos(mutation.getPos());
				seq+=orig.substring(lastPos, curPos)+mutation.getNewNuc();
				lastPos=curPos+1;
			}
			seq+=orig.substring(lastPos);
			return seq;
		}else{
			return orig;
		}
		
	}
	
	public String getMutatedSeq(String chr,int pos, char newNuc){
		if(onChr(chr)&&containsPos(pos)){
			return getMutatedSeq(chrPosToReadPos(pos), newNuc);
		}else{
			return getSeq();
		}
	}
	
	private String getMutatedSeq(int readPos,char newNuc){
		String seq=getSeq();
		return seq.substring(0, readPos)+newNuc+seq.substring(readPos+1);
	}
	
	private int chrPosToReadPos(int chrPos){
		return chrPos-this.pos;
	}
	
	public String toString(){
		return chr+"\t"+pos+"\t"+alignment+"\t"+readID+"\t"+strand+"\t"+mismatches+"\t"+hits+"\t"+read_length+"\t"+offset+"\t"+pe_flag+"\t"+quality;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public String getAlignment() {
		return alignment;
	}

	public String getReadID() {
		return readID;
	}

	public String getStrand() {
		return strand;
	}

	public int getMismatches() {
		return mismatches;
	}

	public int getHits() {
		return hits;
	}

	public int getRead_length() {
		return read_length;
	}

	public int getOffset() {
		return offset;
	}

	public String getPe_flag() {
		return pe_flag;
	}

	public String getQuality() {
		return quality;
	}
	
	
}
