package tools.fastq;

import java.io.BufferedReader;
import java.io.IOException;

public class fastqParser {

	private BufferedReader reader;
	private boolean hasmore= true;
	private String s=null;
	private String currHit= null;
	private String nextHit=null;
	private String currSeq="";
	private String nextSeq=null;
	private String currQual=null;
	private String nextQual=null;
	private String prefix=null;
	
	/**
	 * 
	 * @param reader
	 * @throws IOException
	 * @throws BioException
	 */
	public fastqParser(BufferedReader reader,String prefix) throws IOException, Exception{
		this.reader= reader;
		this.prefix= prefix;
		for(s= reader.readLine();s!=null;s=reader.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='@')
					break;
		readNext();
	}
	public FastqSeq next() throws IOException, Exception{
		String header=nextHit();
		return new FastqSeq(header,getSeq(),getQual());
	}
	public String nextHit() throws IOException, Exception{
		if (!hasmore) {
			throw new Exception("No more hits 1");
		}
		currHit= nextHit;
		currSeq=nextSeq;
		currQual=nextQual;
		readNext();
		return currHit;
	}
	public String getSeq()throws IOException,Exception{
		/*if (currSeq.length()<1) {
			throw new Exception("No more hits 2");
		}*/
		return currSeq;
	}
	public String getQual(){
		return currQual;
	}
//	public String getHit()throws IOException,Exception{
//		/*if (!hasmore) {
//			throw new Exception("No more hits 3");
//		}*/
//		return currHit;
//	}
	public boolean hasNext() {
		return hasmore;
	}
	private void readNext() throws IOException, Exception{
		//assumes that both the sequence and the corresponding quality values only occupy one line each.
		nextSeq="";
		nextQual="";
		if (s==null) {
			hasmore= false;
		} else {
			nextHit=s;
			nextSeq=reader.readLine();
//			s=reader.readLine();
			for(s= reader.readLine();s!=null;s=reader.readLine())
				if(s.length()>0){
					if(s.charAt(0)=='+')
						break;
					else{
						nextSeq+=s;
					}
				}
			if(!nextHit.substring(1).equals(s.substring(1))&&!(s.equals("+"))){
				System.err.println("Strange file:\n"+nextHit+" is not equal to "+s);
			}
			for(s= reader.readLine();s!=null;s=reader.readLine())
				if(s.length()>0){
					if(s.startsWith("@")&&nextQual.length()>=nextSeq.length())
						if(s.startsWith("@"+prefix)){
							break;
						}else{
							System.err.println("quality line:\n"+s);
							nextQual+=s;
						}
					else{
						nextQual+=s;
					}
				}
//			nextQual=reader.readLine();
//			s=reader.readLine();
		}
	}
}
