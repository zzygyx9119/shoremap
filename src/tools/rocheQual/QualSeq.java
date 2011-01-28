package tools.rocheQual;

import java.util.ArrayList;


public class QualSeq{

	private String header,qname;
	private ArrayList<String> seq;

	public QualSeq(String header, String seq){
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
		return this.getSeq(0,this.length());
	}
	
	public String getSeq(int start, int end){
		String s="";
		if(seq.size()>start){
			s+=seq.get(start);
			for(int i=start+1;i<end;i++){
				s+=" "+seq.get(i);
			}
		}
		return s;
	}
	
	public int length(){
		return seq.size();
	}

	public void setSeq(String seq) {
		this.seq= new ArrayList<String>();
		for(String s : seq.trim().split(" +")){
			this.seq.add(s);
		}
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
}
