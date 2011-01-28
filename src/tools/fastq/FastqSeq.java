package tools.fastq;

import java.io.Serializable;
import java.util.ArrayList;

import tools.fasta.FastaSeq;
import tools.rocheQual.RocheQualSeq;

public class FastqSeq implements Serializable{

	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String header,seq,quality;

	public FastqSeq(String header, String seq, String quality) {
		this.header = header;
		this.seq = seq;
		this.quality=quality;
	}
	
	public FastqSeq(FastqSeq fqs){
		this(fqs.getHeader(),fqs.getSeq(),fqs.getQuality());
	}
	
	public FastqSeq(FastaSeq fs){
		this(fs.getHeader(),fs.getSeq(),fs.getGoodSangerQualSeq());
	}

	public String getHeader() {
		return header;
	}
	public String getQname() {
		return header.substring(1).split(" ")[0];
	}

	public void setHeader(String header) {
		this.header = header;
	}

	public String getSeq() {
		return seq;
	}
	
	public String getSeq(int start,int stop){
		return seq.substring(start, stop);
	}
	
	public int length(){
		return seq.length();
	}

	public void setSeq(String seq) {
		this.seq = seq;
	}
	
	public String getQuality() {
		return quality;
	}

	public void setQuality(String quality) {
		this.quality = quality;
	}
	
	public void correct(ArrayList<FastqCorrection> corrections)throws Exception{
		String s="",q="";
		int lastpos= 0;
		for (FastqCorrection correction : corrections) {
			if(this.getSeq().charAt(correction.getPos())!=correction.getFrom()){
				throw new Exception("Trying to apply multiple changes and cannot correct \n"+this.toString()+"\nwith correction:\n"+correction.toString());
			}
			s+=getSeq().substring(lastpos, correction.getPos())+correction.getTo();
			q+=getQuality().substring(lastpos, correction.getPos())+((char)correction.getQual());
			lastpos=correction.getPos()+1;
		}
		this.setSeq(s+getSeq().substring(lastpos));
		this.setQuality(q+getQuality().substring(lastpos));
		
	}
	
	public void correct(FastqCorrection correction)throws Exception{
		if(this.getSeq().charAt(correction.getPos())!=correction.getFrom()){
			throw new Exception("cannot correct \n"+this.toString()+"\nwith correction:\n"+correction.toString());
		}
		this.setSeq(this.getSeq(0, correction.getPos())+correction.getTo()+this.getSeq().substring(correction.getPos()+1));
		this.setQuality(this.getQuality().substring(0, correction.getPos())+((char)correction.getQual())+this.getQuality().substring(correction.getPos()+1));
	}

	public String toString(){
		return header+"\n"+seq+"\n+"+header.substring(1)+"\n"+quality;
	}
	
	public FastaSeq toFastaSeq(){
		return new FastaSeq(">"+this.getHeader().substring(1),this.getSeq());
	}
	
	public RocheQualSeq toPhredQualSeq(){
		return this.toQualSeq(64);
	}
	
//	private String changeQualBase(int origBase, int newBase){
//		String s="";
//		int n;
//		final int diff =origBase-newBase;
//		final String curQuality=this.getQuality();
//		for(int i=0;i<curQuality.length();i++){
//			n=curQuality.charAt(i)-diff;
//			if(n<newBase){
//				System.out.println("Negative quality values when changing base value in: "+ this.getQname());
//				s+=(char) newBase;
//			}else{
//				s+=(char)n;
//			}
//		}
//		return s;
//	}
	
	private String toSangerQual(int base){
		String s="";
		int n;
		for (String number : this.getQuality().split(" +")) {
			if(number.length()>0){
				n=Integer.parseInt(number);
				if(n>=0){
					s+=(char) (base+n);
				}else{
					System.err.println("negative quality values in: "+this.getQname());
					s+=(char) (base);
				}
			}
		}
		return s;
	}
	
	public void convertNumberToPhredQual(){
		setQuality(toSangerQual(64));
	}
	
	private RocheQualSeq toQualSeq(int base){
		String seq="";
		if(this.length()>0){
			seq+=quality.charAt(0)-base;
			for(int i=1;i<this.length();i++){
				seq+=" "+(quality.charAt(i)-base);
			}
		}
		return new RocheQualSeq(">"+header.substring(1), seq);		
	}
	
	private String reverseQual(){
		String revQual="";
		for(int i=0;i<quality.length();i++){
			revQual=quality.charAt(i)+revQual;
		}
		return revQual;
	}
	
	private String reverseComplementSeq(){
		String revCompSeq="";
		for(int i=0;i<seq.length();i++){
			revCompSeq= this.complement(seq.charAt(i))+revCompSeq;
		}
		return revCompSeq;
	}
	
	public FastqSeq reverseComplement(){
		return new FastqSeq(this.header, this.reverseComplementSeq(), this.reverseQual());
	}
	
	public void reverseComplementThis(){
		this.setSeq(this.reverseComplementSeq());
		this.setQuality(this.reverseQual());
	}
	
	private char complement(char c){
		switch (c) {
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'X':
			return 'X';
		case 'N':
			return 'N';
		default:
			System.err.println("unknown nucleotide: "+c+", will return "+c);
			return c;
		}
	}
	
	public String toString(int start,int stop)throws Exception{
		return this.toString(start, stop, header.substring(1));
	}
	
	public String toString(int start,int stop,String newHeader)throws Exception{
		if(start<0||stop>seq.length()){
			throw new Exception("Limits for "+header.substring(1)+", with a length of "+seq.length()+", is out of bounds ("+start+", "+stop+")");
		}
		
		return "@"+newHeader+"\n"+seq.substring(start, stop)+"\n+"+newHeader+"\n"+quality.substring(start,stop);
	}
}
