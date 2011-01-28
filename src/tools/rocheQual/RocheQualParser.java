package tools.rocheQual;

import java.io.BufferedReader;
import java.io.IOException;

public class RocheQualParser {

	private BufferedReader reader;
	private boolean hasmore= true;
	String s=null;
	private String currHit= null;
	private String nextHit=null;
	private String currSeq="";
	private String nextSeq=null;
	
	/**
	 * 
	 * @param reader
	 * @throws IOException
	 * @throws BioException
	 */
	public RocheQualParser(BufferedReader reader) throws IOException, Exception{
		this.reader= reader;
		for(s= reader.readLine();s!=null;s=reader.readLine())
			if(s.length()>0)
				if(s.charAt(0)=='>')
					break;
		readNext();
	}
	public RocheQualSeq next() throws IOException, Exception{
		String header=nextHit();
		return new RocheQualSeq(header,getSeq());
	}
	public String nextHit() throws IOException, Exception{
		if (!hasmore) {
			throw new Exception("No more hits 1");
		}
		currHit= nextHit;
		currSeq=nextSeq;
		readNext();
		return currHit;
	}
	public String getSeqOrig()throws IOException,Exception{
		/*if (currSeq.length()<1) {
			throw new Exception("No more hits 2");
		}*/
		return currSeq;
	}
	public String getSeq()throws IOException,Exception{
		/*if (currSeq.length()<1) {
			throw new Exception("No more hits 2");
		}*/
		return currSeq.replace("\n", "");
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
		nextSeq="";
		if (s==null) {
			hasmore= false;
		} else {
			nextHit=s;
			for(s= reader.readLine();s!=null;s=reader.readLine())
				if(s.length()>0){
					if(s.charAt(0)=='>')
						break;
					else{
						nextSeq+=" "+s;
					}
				}
		}
	}
}
