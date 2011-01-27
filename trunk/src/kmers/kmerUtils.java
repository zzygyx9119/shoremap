package kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class kmerUtils {

	private static HashMap<String, String> methodsHelp= new HashMap<String, String>();
	private final static String sep="\t";
	private final static String[] nucleotides= new String[]{"A","C","G","T","N"};
	
	public static void main(String[] args)throws Exception{
		//load helpMap
		methodsHelp.put("sortKmers", "sortKmers - sorts the kmers\n\targs = <kmer file> <kmer position> <bin size> <tmpPrefix>\n");
		methodsHelp.put("sortKmersLink", "sortKmersLink - sorts the kmers and adds a flag if it is linked to the previous read\n\targs = <kmer file> <kmer position> <max jump>\n");
		methodsHelp.put("generateKmers", "generateKmers - generates all kmers of length n\n\targs = <n>\n");
		methodsHelp.put("generateSeeds", "generateSeeds - generates seeds from a kmer file\n\targs = <kmerFile> <kmer position> <min seed size> <outPrefix>\n");
		
		//check which method to run
		if(args.length>0){
			if(args[0].equals("sortKmers")&&args.length==5){
				sortKmers(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),args[4]);
			}else if(args[0].equals("sortKmersLink")&&args.length==4){
				sortKmersLink(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("generateKmers")&&args.length==2){
				generateKmers(Integer.parseInt(args[1]),"");
			}else if(args[0].equals("generateSeeds")&&args.length==5){
				generateSeeds(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),args[4]);
			}else if(args[0].equals("method2")&&args.length>1){
				
			}else{
				System.err.println(printHelp(args[0]));
			}
		}else{
			System.err.println(printHelp());
		}
	}
	
	private static String printHelp(){
		String help="";
		for (String s : methodsHelp.values()) {
			help+=s;
		}
		
		return help;
	}
	
	private static String printHelp(String method){
		if(methodsHelp.containsKey(method)){
			return methodsHelp.get(method);
		}
		return printHelp();
	}
	
	private static void generateKmers(int n,String s)throws Exception{
		if(n==0){
			System.out.println(s);
		}else{
			for(int i=0;i<4;i++){
				generateKmers(n-1, nucleotides[i]+s);
			}
		}
	}
	
	private static void generateSeeds(String kmerFile,int pos,int minSeedSize,String outPrefix)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		ArrayList<kmerData> kmers= new ArrayList<kmerData>();
		HashSet<String> suffixes= new HashSet<String>(1000000);
		kmerData kmer;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			kmer=new kmerData(pos, s);
			kmers.add(kmer);
			suffixes.add(kmer.suffix());
		}
		Collections.sort(kmers);
		//Find seeds
//		seedData seed= new seedData(kmers.get(0));
		int seedCount=0;
//		boolean exist;
		ArrayList<seedData> seeds= new ArrayList<seedData>();
//		HashSet<String> usedKmers= new HashSet<String>(1000000);
		BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+".fa"));
		BufferedWriter outLong= new BufferedWriter(new FileWriter(outPrefix+"_long.txt"));
		
		boolean found=true;
		while(found){
			System.err.print("\t"+kmers.size()+" seeding"+"\r");
			found=false;
			//find unique seeds
			seeds=new ArrayList<seedData>();
			for(int i=0;i<kmers.size();i++){
				kmer=kmers.get(i);
				if(!suffixes.contains(kmer.prefix())){
					seeds.add(new seedData(kmer));
				}
			}
			System.err.print("\t"+kmers.size()+seeds.size()+"      "+"\r");
			found=seeds.size()>0;
			//extend seeds
			HashSet<kmerData> usedKmers=new HashSet<kmerData>();
			for (seedData initialSeed : seeds) {
				for (seedData finalSeed : generateSeeds_extend(initialSeed, kmers)) {
					usedKmers.addAll(finalSeed.getKmers());
					if(finalSeed.size()>=minSeedSize){
						seedCount++;
						out.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.consensus()+"\n");
						outLong.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.toString()+"\n");
					}
				}
			}
			//clean
			kmers.removeAll(usedKmers);
			for (kmerData kmerData : usedKmers) {
				suffixes.remove(kmerData.suffix());
			}
		}
		while(kmers.size()>0){
			System.err.println("\t"+kmers.size()+"                \r");
			//extract cyclic seeds
			for(seedData finalSeed : generateSeeds_extend(new seedData(kmers.get(0)), kmers)){
				kmers.removeAll(finalSeed.getKmers());
				if(finalSeed.size()>=minSeedSize){
					seedCount++;
					out.write(">"+seedCount+" "+finalSeed.size()+" cyclic\n"+finalSeed.consensus()+"\n");
					outLong.write(">"+seedCount+" "+finalSeed.size()+" cyclic\n"+finalSeed.toString()+"\n");
				}
			}
		}
//		
//		
//		for(int i=1;i<kmers.size();i++){
//			if(i%100==0){
//				System.err.print("\t"+i+" "+seeds.size()+"\r");
//				out.flush();
//				outLong.flush();
//			}
//			kmer=kmers.get(i);
//			if(!usedKmers.contains(kmer.getKmer())){
//				if(!seed.addLast(kmer)){
//					exist=false;
//					for(int j=0;j<seeds.size()&&!exist;j++){
//						if(seeds.get(j).contains(seed)){
//							exist=true;
//						}
//					}
//					for (seedData finalSeed : generateSeeds_extend(seed,kmers)){
//						for (kmerData kd : finalSeed.getKmers()) {
//							usedKmers.add(kd.getKmer());
//						}
//						if(finalSeed.size()>=minSeedSize){
//							seedCount++;
//							out.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.consensus()+"\n");
//							outLong.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.toString()+"\n");
//						}
//					}
//					if(!exist){
//						seeds.addAll(generateSeeds_extend(seed,kmers));
//					}
//					seed=new seedData(kmers.get(i));
//				}
//			}
//		}
//		for (seedData finalSeed : seeds) {
//			if(finalSeed.size()>=minSeedSize){
//				seedCount++;
//				out.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.consensus()+"\n");
//				outLong.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.toString()+"\n");
//			}
//		}
		out.close();
		outLong.close();
//		return seeds;
	}
	
	private static ArrayList<seedData> generateSeeds_extend(seedData seed, ArrayList<kmerData> kmers){
		seedData newSeed= new seedData(seed);
		ArrayList<seedData> seeds= new ArrayList<seedData>();
		boolean fin=true;
		for(int i=0;i<kmers.size();i++){
			if(newSeed.addLast(kmers.get(i))){
				seeds.addAll(generateSeeds_extend(newSeed,kmers));
				newSeed= new seedData(seed);
				fin=false;
			}
		}
		if(fin){
			seeds.add(seed);
		}
		return seeds;
	}
	
	
	private static void sortKmersLink(String kmerFile, int pos, int maxJump)throws Exception{
//		System.err.println("kmerFile:"+kmerFile+"\npos: "+pos+"\nmaxJump: "+maxJump);
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		ArrayList<kmerData> kmers= new ArrayList<kmerData>();
//		System.err.println("Reading "+kmerFile);
		int line=0;
		for(String s=in.readLine();s!=null;s=in.readLine()){
//			if(line%10000==0){
//				System.out.print("\t"+line+"\r");
//			}
			kmers.add(new kmerData(pos, s));
			++line;
		}
//		System.err.println(line+"\nsorting...");
		Collections.sort(kmers);
//		System.err.println("printing and idenifying links");
		String lastKmer="XXXXXXXXXXX";
		while(lastKmer.length()<maxJump+1);{
			lastKmer+="X";
		}
		boolean linked;
		for (kmerData kmerData : kmers) {
			System.out.print(kmerData);
			linked=false;
			for(int i=1;i<=maxJump&&!linked;i++){
				linked=kmerData.getKmer().startsWith(lastKmer.substring(i));
			}
			if(linked){
				System.out.println(sep+"1");
			}else{
				System.out.println(sep+"0");
			}
			lastKmer=kmerData.getKmer();
		}
	}
	
	private static void sortKmers(String kmerFile,int pos, int binSize,String tmpPrefix)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		ArrayList<kmerData> kmers= new ArrayList<kmerData>();
		for(String s=in.readLine();s!=null;s=in.readLine()){
			kmers.add(new kmerData(pos, s));
			if(kmers.size()>=binSize){
				Collections.sort(kmers);
				//push to disk
				
			}
		}
		Collections.sort(kmers);
		//do mergeSort...
		
		for (kmerData kmerData : kmers) {
			System.out.println(kmerData);
		}
	}
}

class seedData{
	
	private ArrayList<kmerData> kmers;
	
	private seedData(){
		kmers= new ArrayList<kmerData>();
	}
	
	public seedData(kmerData first){
		this();
		kmers.add(first);
	}
	
	public seedData(seedData old){
		this();
		for (kmerData kd : old.kmers) {
			this.kmers.add(kd);
		}
	}
	
	public ArrayList<kmerData> getKmers(){
		return kmers;
	}
	
	public String consensus(){
		//presumes that the kmers are of the same length
		String seq=kmers.get(0).getKmer();
		int lastPos=seq.length()-1;
		for(int i=1;i<kmers.size();i++){
			seq+=""+kmers.get(i).getKmer().charAt(lastPos);
		}
		return seq;
	}
	
	public boolean addLast(kmerData ext){
		if(extendEnd(ext)&&!contains(ext)){
			kmers.add(ext);
			return true;
		}
		return false;
	}
	
	public boolean addFirst(kmerData ext){
		if(extendStart(ext)&&!contains(ext)){
			kmers.add(0, ext);
			return true;
		}
		return false;
	}
	
	public boolean extendEnd(kmerData ext){
		return kmers.get(kmers.size()-1).endsWith(ext);
	}
	
	public boolean extendStart(kmerData ext){
		return kmers.get(0).startsWith(ext);
	}
	
	public boolean contains(kmerData kd){
		return kmers.contains(kd);
	}
	
	public boolean contains(seedData seed){
		return this.kmers.containsAll(seed.kmers);
	}
	
	public int size(){
		return kmers.size();
	}
	
	public String toString(){
		String out=consensus();
		for (kmerData kmer : kmers) {
			out+="\n"+kmer.toString();
		}
		return out;
	}
}

class kmerData implements Comparable<kmerData>,Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int pos;
	private int permutations;
	private String rep;
	private String[] data;
	
	public kmerData(int pos, String s) throws Exception{
		this(pos, s.split("\t"));
	}
	
	public kmerData(int pos, String[] data) throws Exception {
		super();
		this.pos = pos;
		this.data = data;
		rep=data[pos];
		permutations=0;
		String tmp=rep;
		for(int i=0;i<rep.length();i++){
			tmp=tmp.substring(1)+tmp.charAt(0);
			if(tmp.compareTo(rep)<=0){
				rep=tmp;
				permutations=0;
			}else{
				permutations++;
			}
			if(data[pos].equals(tmp)){
				break;
			}
		}
	}
	
	public String suffix(int n){
		return getKmer().substring(n);
	}
	
	public String suffix(){
		return suffix(1);
	}
	
	public String prefix(int n){
		return getKmer().substring(0, getKmer().length()-n);
	}
	
	public String prefix(){
		return prefix(1);
	}
	
	public boolean startsWith(String prefix){
		return getKmer().startsWith(prefix);
	}
	
	public boolean startsWith(kmerData prefix,int n){
		return startsWith(prefix.suffix(n));
//		return startsWith(prefix.getKmer().substring(n));
	}
	
	public boolean startsWith(kmerData prefix){
		return startsWith(prefix,1);
	}
	
	public boolean endsWith(String suffix){
		return getKmer().endsWith(suffix);
	}
	
	public boolean endsWith(kmerData suffix, int n){
		return endsWith(suffix.prefix(n));
		//return endsWith(suffix.getKmer().substring(0, suffix.getKmer().length()-n));
	}
	
	public boolean endsWith(kmerData suffix){
		return endsWith(suffix,1);
	}
	
	public String getKmer(){
		return data[pos];
	}

	@Override
	public int compareTo(kmerData o) {
		if(this.rep.equals(o.rep))
			return this.permutations-o.permutations;
		return this.rep.compareTo(o.rep);
	}
	
	public String toString(){
		return this.toString("\t");
	}
	
	public String toString(String sep){
		String s= data[0];
		for(int i=1;i<data.length;i++){
			s+=sep+data[i];
		}
		return s;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(data);
		result = prime * result + pos;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		kmerData other = (kmerData) obj;
		if (!Arrays.equals(data, other.data))
			return false;
		if (pos != other.pos)
			return false;
		return true;
	}
	
	
}
