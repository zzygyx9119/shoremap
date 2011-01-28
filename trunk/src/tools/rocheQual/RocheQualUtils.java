package tools.rocheQual;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;


public class RocheQualUtils {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
		if(args.length>0){
			if(args[0].equals("tag")&&args.length>2){
				tag(args[1],args[2]);
			}else if(args[0].equals("antisub")&&args.length>2){
				sub(args[1],args[2],false);
			}else if(args[0].equals("sub")&&args.length>2){
				sub(args[1],args[2],true);
			}else if(args[0].equals("gmToStrainGFF")&&args.length>3){
				//gmToStrainGFF(args[1], args[2], args[3]);
			}else{
				System.err.println(printHelp());
				System.exit(616);
			}
		}else{
			System.err.println(printHelp());
			System.exit(616);
		}
	}

	private static String printHelp(){
		String help="Usage: qualUtils <cmd> <input>\n";
		help+="where <cmd> is:\n";
		help+="tag - tags the qual header: >TAG_qname\n";
		help+="\t<input> = <qual file> <TAG>\n";
		help+="sub - extracts the sequences in the list\n";
		help+="\t<input> = <qual file> <giList file>\n";
		help+="antisub - extracts the sequences NOT in the list\n";
		help+="\t<input> = <qual file> <giList file>\n";
				
		
		return help;
	}
	
	public static void sub(String qualFile,String giListFile,boolean sub)throws Exception{
		//ArrayList<String> giList= new ArrayList<String>();
		HashSet<String> giList= new HashSet<String>();
		
		BufferedReader in= new BufferedReader(new FileReader(giListFile));
		for(String s=in.readLine();s!=null;s=in.readLine()){
			giList.add(s);
		}
		RocheQualParser qp =new RocheQualParser(new BufferedReader(new FileReader(qualFile)));
		RocheQualSeq qs;
		for(;qp.hasNext();){
			qs=qp.next();
			if(giList.contains(qs.getQname())==sub){
				System.out.println(qs);
			}
		}
	}
	
	public static void tag(String qualFile,String tag)throws Exception{
		RocheQualParser qp =new RocheQualParser(new BufferedReader(new FileReader(qualFile)));
		for(;qp.hasNext();){
			String s=qp.nextHit();
			System.out.println(">"+tag+"_"+s.substring(1)+"\n"+qp.getSeq());
		}
	}
}
