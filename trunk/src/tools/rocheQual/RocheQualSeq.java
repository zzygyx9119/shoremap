package tools.rocheQual;

import tools.fasta.FastaSeq;




public class RocheQualSeq{

	// this is the roche qual seq with numbers instead of ascii chars
	
	private String header,qname;
	//private ArrayList<String> seq;
	private String seq;

	public RocheQualSeq(String header, String seq){
		this.header = header;
		this.updateQname();
		this.setSeq(seq);
	}
	
	

	public String getHeader() {
		return header;
	}

	public void setHeader(String header) {
		this.header = header;
		this.updateQname();
	}

	public String getSeq() {
		return this.seq;
	}
	
	public String getSeq(int start, int end){
		String s="";
		if(this.length()>start){
			final String[] l=s.trim().split(" +");
			s+=l[start];
			for(int i=start+1;i<end;i++){
				s+=" "+l[i];
			}
		}
		return s;
	}
	
	public int[] getIntSeq(){
		String[] l=seq.trim().split(" +");
		int[] quals=new int[l.length];
		for(int i=0;i<quals.length;i++){
			quals[i]=Integer.parseInt(l[i]);
		}
		return quals;
	}
	
	public int length(){
		return seq.trim().split(" +").length;
	}

	public void setSeq(String seq) {
		this.seq= seq;
	}
	
	public String getQname(){
		return qname;
	}
	
	public String toString(){
		return header+"\n"+this.getSeq();
	}
	
	public String toString(int start,int stop)throws Exception{
		if(start<0||stop>this.length()){
			throw new Exception("Limits for "+qname+", with a length of "+this.length()+", is out of bounds ("+start+", "+stop+")");
		}
		
		return header+"\n"+this.getSeq(start, stop);
	}
	
	private void updateQname(){
		qname=header.substring(1).split(" ")[0];
	}
	
	public FastaSeq	toSangerQual(){
		return new FastaSeq(this.getHeader(), this.toSangerQual(33));
	}
	
	private String toSangerQual(int base){
		String s="";
		int n;
		for (String number : this.getSeq().split(" +")) {
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
}
