package tools.fasta;

public class FastaSeq {

	private String header,seq,qname;

	public FastaSeq(String header, String seq) {
		this.header = header;
		this.updateQname();
		this.seq = seq;
	}
	
	public String getGoodSangerQualSeq(){
		return getSangerQualSeq(30);
	}
	
	private String getSangerQualSeq(int quality){
		return getQualSeq(quality,33);
	}
	
	private String getQualSeq(int quality, int base){
		String qual="";
		for(int i=0;i<this.length();i++){
			qual+=(char) (base+quality);
		}
		return qual;
	}

	public String getHeader() {
		return header;
	}

	public void setHeader(String header) {
		this.header = header;
		this.updateQname();
	}

	public String getSeq() {
		return seq;
	}
	
	public String getSeq(int start,int end){
		return seq.substring(start, end);
	}
	
	public int length(){
		return seq.length();
	}

	public void setSeq(String seq) {
		this.seq = seq;
	}
	
	public String getQname(){
		return qname;
	}
	
	public String toString(){
		return header+"\n"+seq;
	}
	
	public String toString(int start,int stop)throws Exception{
		if(start<0||stop>seq.length()){
			throw new Exception("Limits for "+qname+", with a length of "+seq.length()+", is out of bounds ("+start+", "+stop+")");
		}
		
		return header+"\n"+this.getSeq(start, stop);
	}
	
	private void updateQname(){
		qname=header.substring(1).split(" ")[0];
	}
	
}
