package tools.reptile;

import java.util.ArrayList;

import tools.fastq.FastqCorrection;

public class ReptileLine {

	final char[] alphabet= new char[] {'A','C','G','T','N'};
	int read;
	ArrayList<FastqCorrection> corrections;
	
	public ReptileLine(String line){
		String[] l=line.split("\t");
		read=Integer.parseInt(l[0]);
		corrections= new ArrayList<FastqCorrection>();
		for(int i=2;i<l.length;i+=4){
			corrections.add(new FastqCorrection(Integer.parseInt(l[i]), Integer.parseInt(l[i+3]), alphabet[Integer.parseInt(l[i+1])], alphabet[Integer.parseInt(l[i+2])]));
		}
	}
	
	public ArrayList<FastqCorrection> getCorrections(){
		return corrections;
	}
	
	public int getReadId(){
		return read;
	}
}
