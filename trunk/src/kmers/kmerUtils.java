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

import tools.fastq.fastqParser;
import tools.fastq.fastqUtils;
import tools.shoreMap.Mutation;
import tools.shoreMap.ShoreMapLine;

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
		methodsHelp.put("mapListMutationKmers", "mapListMutationKmers - takes a map.list and a list of mutations (Chr\\tposition\\tnewNuc) and prints two files, one with the kmers of size n to remove and one with the ones to add\n\targs = <map.list> <mutation file> <kmer size n> <outPrefix>\n");
		methodsHelp.put("kmerize", "kmerize - prints all kmers of size n in each sequence\n\targs = <fastqFile> <n>\n");
		methodsHelp.put("kmerizeEnd", "kmerizeEnd - prints the kmers in the end of a read missing when reducing from n to m\n\targs = <fastqFile> <n> <m>\n");
		
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
			}else if(args[0].equals("mapListMutationKmers")&&args.length==5){
				mapListMutationKmers(args[1],args[2],Integer.parseInt(args[3]),args[4]);
			}else if(args[0].equals("kmerize")&&args.length==3){
				kmerize(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerize")&&args.length==3){
				kmerize(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerizeEnd")&&args.length==4){
				kmerizeEnd(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
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
	
	protected static String reverseComplementSeq(String seq){
		String revCompSeq="";
		for(int i=0;i<seq.length();i++){
			revCompSeq= complement(seq.charAt(i))+revCompSeq;
		}
		return revCompSeq;
	}
	
	protected static String reverse(String seq){
		String rev="";
		for(int i=0;i<seq.length();i++){
			rev=seq.charAt(i)+rev;
		}
		return rev;
	}
	
	protected static String complement(String seq){
		String comp="";
		for(int i=0;i<seq.length();i++){
			comp=comp+complement(seq.charAt(i));
		}
		return comp;
	}
	
	private static char complement(char c){
		switch (c) {
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		case 'N':
			return 'N';
		default:
			System.err.println("unknown nucleotide: "+c+", will return "+c);
			return c;
		}
	}
	
	protected static String kmerToUse(String kmer){
		String revKmer=reverseComplementSeq(kmer);
		if(kmer.compareTo(revKmer)<0){
			return kmer;
		}else{
			return revKmer;
		}
	}
	
	private static void printKmer(String kmer){
		System.out.println(kmerToUse(kmer));
	}
	
	public static void kmerize(String fastqFile, int n)throws Exception{
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),fastqUtils.getPrefix(fastqFile));
		String seq;
		for(;fqp.hasNext();){
			seq=fqp.next().getSeq();
			for(int i=0,j=n;j<=seq.length();i++,j++){
				printKmer(seq.substring(i, j));
			}
		}
	}
	
	public static void kmerizeEnd(String fastqFile, int n, int m)throws Exception{
		if(m>=n){
			throw new Exception("m ("+m+") must be smaller than n ("+n+")");
		}
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),fastqUtils.getPrefix(fastqFile));
		for(;fqp.hasNext();){
			String seq= fqp.next().getSeq();
			seq=seq.substring(seq.length()-n+1);
			for(int i=0,j=m;j<=seq.length();i++,j++){
				printKmer(seq.substring(i, j));
			}
		}
	}
	
	private static void mapListMutationKmers(String mapListFile,String mutationFile,int n, String outPrefix) throws Exception{
		BufferedReader mutationReader = new BufferedReader(new FileReader(mutationFile));
		HashMap<String, ArrayList<Mutation>> mutations= new HashMap<String, ArrayList<Mutation>>();
		Mutation mutation;
		for(String s=mutationReader.readLine();s!=null;s=mutationReader.readLine()){
			String l[]=s.split("\t");
			mutation=new Mutation(l[0], Integer.parseInt(l[1]), l[2].charAt(0));
			if(!mutations.containsKey(mutation.getChr())){
				mutations.put(mutation.getChr(), new ArrayList<Mutation>());
			}
			mutations.get(mutation.getChr()).add(mutation);
		}
		mutationReader.close();
		BufferedReader mapListReader = new BufferedReader(new FileReader(mapListFile));
		ShoreMapLine sml;
		BufferedWriter origWriter= new BufferedWriter(new FileWriter(outPrefix+"_orig.kmers"));
		BufferedWriter mutatedWriter= new BufferedWriter(new FileWriter(outPrefix+"_mutated.kmers"));
		for(String s= mapListReader.readLine();s!=null;s=mapListReader.readLine()){
			sml= new ShoreMapLine(s);
			if(mutations.containsKey(sml.getChr())){
				//check for mutation... this has to be rewritten if there are a lot of mutations or an unfiltered map.list
				String orig=sml.getSeq();
				String mutated=sml.getMutatedSeq(mutations.get(sml.getChr()));
				if(!orig.equals(mutated)){
					//kmerize
					for(int i=0,j=n;j<orig.length();i++,j++){
						origWriter.write(kmerToUse(orig.substring(i,j))+"\n");
						mutatedWriter.write(kmerToUse(mutated.substring(i, j))+"\n");
					}
				}
			}
		}
		mapListReader.close();
		origWriter.close();
		mutatedWriter.close();
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
			suffixes.add(kmerToUse(kmer.suffix()));
		}
//		Collections.sort(kmers);
		//Find seeds
//		seedData seed= new seedData(kmers.get(0));
		int seedCount=0;
//		boolean exist;
		ArrayList<seedData> seeds= new ArrayList<seedData>();
//		HashSet<String> usedKmers= new HashSet<String>(1000000);
		BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+".fa"));
		BufferedWriter outLong= new BufferedWriter(new FileWriter(outPrefix+"_long.txt"));
		
		boolean found=true;
		System.err.println("Unique seeds...");
		while(found){
			System.err.print("    "+kmers.size()+" seeding"+"\r");
			found=false;
			//find unique seeds
			seeds=new ArrayList<seedData>();
			for(int i=0;i<kmers.size();i++){
				kmer=kmers.get(i);
				if(!suffixes.contains(kmerToUse(kmer.prefix()))){
					seeds.add(new seedData(kmer));
				}
			}
			System.err.print("    "+kmers.size()+" "+seeds.size()+"      "+"\r");
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
				suffixes.remove(kmerToUse(kmerData.suffix()));
			}
		}
		System.err.println("Repetitive structures...");
		while(kmers.size()>0){
			System.err.print("    "+kmers.size()+"                \r");
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
		System.err.println("Done!");
	}
	
	private static ArrayList<seedData> generateSeeds_extend(seedData seed, ArrayList<kmerData> kmers){
		seedData newSeed= new seedData(seed);
		ArrayList<seedData> seeds= new ArrayList<seedData>();
		boolean fin=true;
		for(int i=0;i<kmers.size();i++){
			if(newSeed.addRight(kmers.get(i))){
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
	
	public boolean addRight(kmerData ext){
		if(extendRight(ext)&&!contains(ext)){
			kmers.add(ext);
			return true;
		}
		return false;
	}
	
	public boolean addLeft(kmerData ext){
		if(extendLeft(ext)&&!contains(ext)){
			kmers.add(0, ext);
			return true;
		}
		return false;
	}
	
	public boolean extendRight(kmerData ext){
		return kmers.get(kmers.size()-1).extendRight(ext);
	}
	
	public boolean extendLeft(kmerData ext){
		return kmers.get(0).extendLeft(ext);
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
		String out=consensus()+"\n"+kmers.get(0);
		String lastSuffix=kmers.get(0).suffix();
		kmerData kmer;
		for(int i=1;i<kmers.size();i++){
			kmer=kmers.get(i);
			if(kmer.getKmer().startsWith(lastSuffix)){
				out+="\n"+kmer.toString();
			}else{
				out+="\n"+kmer.toStringRev();
			}
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
	private String kmerSeqComp;
	private String[] data;
	
	public kmerData(int pos, String s) throws Exception{
		this(pos, s.split("\t"));
	}
	
	public kmerData(int pos, String[] data) throws Exception {
		super();
		this.pos = pos;
		this.data = data;
		kmerSeqComp=kmerUtils.complement(data[pos]);
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
	
//	public String suffixRev(){
//		return suffixRev(1);
//	}
//	
//	public String suffixRev(int n){
//		return getKmerRev().substring(n);
//	}
//	
//	public String prefixRev(int n){
//		return getKmerRev().substring(0, getKmer().length()-n);
//	}
//	
//	public String prefixRev(){
//		return prefix(1);
//	}
	
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
	
	private boolean extendLeft(String prefix){
		return getKmer().startsWith(prefix)||getKmerComp().startsWith(prefix);
	}
	
	public boolean extendLeft(kmerData prefix,int n){
		return extendLeft(prefix.suffix(n));
//		return startsWith(prefix.getKmer().substring(n));
	}
	
	public boolean extendLeft(kmerData prefix){
		return extendLeft(prefix,1);
	}
	
	private boolean extendRight(String suffix){
		return getKmer().endsWith(suffix)||getKmerComp().endsWith(suffix);
	}
	
	public boolean extendRight(kmerData suffix, int n){
		return extendRight(suffix.prefix(n));
		//return endsWith(suffix.getKmer().substring(0, suffix.getKmer().length()-n));
	}
	
	public boolean extendRight(kmerData suffix){
		return extendRight(suffix,1);
	}
	
	public String getKmer(){
		return data[pos];
	}
	
	private String getKmerComp(){
		return kmerSeqComp;
	}

	public int compareTo(kmerData o) {
		if(this.rep.equals(o.rep))
			return this.permutations-o.permutations;
		return this.rep.compareTo(o.rep);
	}
	
	public String toString(){
		return this.toString("\t");
	}
	
	public String toStringRev(){
		return toStringRev("\t");
	}
	
	public String toStringRev(String sep){
		String s="";
		if(pos==0){
			s=getKmerComp();
		}else{
			s=data[0];
		}
		for(int i=1;i<data.length;i++){
			if(i==pos){
				s=sep+getKmerComp();
			}else{
				s=sep+data[i];
			}
		}
		return s;
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
