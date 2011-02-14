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

import tools.fastq.FastqSeq;
import tools.fastq.fastqParser;
import tools.kmer.KmerSet_binary;
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
		methodsHelp.put("reduce", "reduce - reduces the kmers at pos in a file to length n\n\targs = <fastqFile> <pos> <n>\n");
		methodsHelp.put("extractGood", "extractGood - Takes a kmerFile with the kmers at pos (column count starts at 0), and extract all fastq seqs that contains at least min of the good kmers\n\targs = <kmerFile> <pos> <fastqFile> <min> <outPrefix>\n");
		
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
			}else if(args[0].equals("reduce")&&args.length==4){
				reduce(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("kmerize")&&args.length==3){
				kmerize(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerizeEnd")&&args.length==4){
				kmerizeEnd(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("extractGood")&&args.length==6){
				extractGood(readKmers(args[1],Integer.parseInt(args[2])),args[3],Integer.parseInt(args[4]),args[5]);
			}else if(args[0].equals("method2")&&args.length>1){
				
			}else{
				System.err.println(printHelp(args[0]));
			}
		}else{
			System.err.println(printHelp());
		}
	}
	
	public static void extractGood(KmerSet_binary kmers,String fastqFile, int min, String outPrefix)throws Exception{
		BufferedWriter fastqout= new BufferedWriter(new FileWriter(outPrefix+"_cleaned.fastq"));
		BufferedWriter countout= new BufferedWriter(new FileWriter(outPrefix+"_count.csv"));
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)), "");
		FastqSeq fqs;
		int kmerSize=kmers.getKmerSize(),count;
		System.err.println("Parsing fastq file...");
		for(int i=0;fqp.hasNext();++i){
			if(i%10000==0){
				System.err.print("    "+i+"           \r");
			}
			fqs=fqp.next();
			count=0;
			for (String kmer : kmerizeString(fqs.getSeq(),kmerSize)) {
				if(kmers.exists(kmer)){
					++count;
				}
			}
			if(count>=min){
				fastqout.write(fqs.toString()+"\n");
				countout.write(fqs.getQname()+"\t"+count+"\n");
			}
		}
		
		fastqout.close();
		countout.close();
	}
	
	private static ArrayList<String> kmerizeString(String seq, int length){
		ArrayList<String> seqs = new ArrayList<String>();
		if(seq.length()>=length){
			String kmer=seq.substring(0, length);
			seqs.add(kmerToUse(kmer));
			for(int i=length;i<seq.length();i++){
				kmer=kmer.substring(1)+seq.charAt(i);
				seqs.add(kmerToUse(kmer));
			}
		}
		return seqs;
	}
	
	private static KmerSet_binary readKmers(String kmerCountFile, int pos)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerCountFile));
		String[] l=in.readLine().split("\t");
		KmerSet_binary kmers= new KmerSet_binary(l[pos]);
		int line=1;
		for(String s=in.readLine();s!=null;s=in.readLine(),line++){
			if(line%10000==0)
				System.err.print("\t"+line+"   "+kmers.size()+"\r");
			l=s.split("\t");
			kmers.addSeq(l[pos]);
		}
		
		return kmers;
	}
	
	public static void reduce(String kmerFile,int pos,int length)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerFile));
		for(String s=in.readLine();s!=null;s=in.readLine()){
			System.out.println((new kmerData(pos, s)).reduce(length));
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
	
//	protected static String kmerToUse(String kmer){
//		String revKmer=reverseComplementSeq(kmer);
//		if(kmer.compareTo(revKmer)<0){
//			return kmer;
//		}else{
//			return revKmer;
//		}
//	}
	
	public static String kmerToUse(String kmer){
		String kmerRev= "";
		int j=kmer.length()-1,result;
		char c,r;
		for(int i=0;i<kmer.length();i++,j--){
			c=kmer.charAt(i);
			r=complement(kmer.charAt(j));
			if((result=c-r)!=0){
				if(result<0){
					return kmer;
				}else{
					kmerRev+=""+r;
					++i;
					break;
				}
			}else{
				kmerRev=r+kmerRev;
			}
		}
		for(;j>=0;j--){
			kmerRev+=complement(kmer.charAt(j))+"";
		}
		return kmerRev;
	}
	
	private static void printKmer(String kmer){
		System.out.println(kmerToUse(kmer));
	}
	
	public static void kmerize(String fastqFile, int n)throws Exception{
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
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
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
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
		HashMap<String,HashSet<kmerData>> kmers= new HashMap<String, HashSet<kmerData>>();
		HashSet<String> suffixes= new HashSet<String>(1000000);
//		HashSet<String> prefixes= new HashSet<String>(1000000);
		
		HashSet<kmerData> uniqueKmers= new HashSet<kmerData>();
		
		kmerData kmer;
		for(String s=in.readLine();s!=null;s=in.readLine()){
			kmer=new kmerData(pos, s);
			if(!kmers.containsKey(kmer.prefix())){
				kmers.put(kmer.prefix(), new HashSet<kmerData>());
			}
			if(!kmers.containsKey(kmer.prefixRevComp())){
				kmers.put(kmer.prefixRevComp(), new HashSet<kmerData>());
			}
			kmers.get(kmer.prefix()).add(kmer);
			kmers.get(kmer.prefixRevComp()).add(kmer);
			suffixes.add(kmer.suffix());
			suffixes.add(kmer.suffixRevComp());
			uniqueKmers.add(kmer);
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
		HashSet<kmerData> usedKmers;
		boolean found=true;
		System.err.println("Unique seeds...");
//		HashSet<String> printedSeeds= new HashSet<String>();
//		while(kmers.size()>0){
//			System.err.print("    "+kmers.size()+"\r");
//
//			HashSet<kmerData> usedKmers=new HashSet<kmerData>();
//			
//			for (seedData finalSeed : generateSeeds_extend(new seedData(kmers.get(0)), kmers)) {
//				usedKmers.addAll(finalSeed.getKmers());
//				if(finalSeed.size()>=minSeedSize&&printedSeeds.add(finalSeed.consensus())){
//					seedCount++;
//					out.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.consensus()+"\n");
//					outLong.write(">"+seedCount+" "+finalSeed.size()+"\n"+finalSeed.toString()+"\n");
//				}
//			}
//			kmers.removeAll(usedKmers);
//		}
		while(found){
			found=false;			
			System.err.print("    "+kmers.size()+"\r");
			//find unique seeds
			seeds=new ArrayList<seedData>();
			for (kmerData kmerData : uniqueKmers) {
				if(!suffixes.contains(kmerData.prefix())){
					seeds.add(new seedData(kmerData));
				}
//				if(!suffixes.contains(kmerData.prefixRevComp())){
//					seeds.add(new seedData(kmerData.revComp()));
//				}
			}
			System.err.print("    "+kmers.size()+" "+seeds.size()+"      "+"\r");
			found=seeds.size()>0;
			//extend seeds
			usedKmers=new HashSet<kmerData>();
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
			uniqueKmers.removeAll(usedKmers);
			
			for (kmerData kmerData : usedKmers) {
				kmers.remove(kmerData.prefix());
				kmers.remove(kmerData.prefixRevComp());
			}
			//Rebuild the suffix list
			suffixes= new HashSet<String>(1000000);
			for (kmerData kmerData : uniqueKmers) {
				suffixes.add(kmerData.suffix());
				suffixes.add(kmerData.suffixRevComp());
			}
		}
		System.err.println("Repetitive structures...");
		while(kmers.size()>0){
//			System.err.print("    "+kmers.size()+"                \r");
			//extract cyclic seeds
			HashSet<kmerData> onePrefix= new HashSet<kmerData>();
			for (HashSet<kmerData> tmp : kmers.values()) {
				onePrefix=new HashSet<kmerData>(tmp);
				break;
			}
			usedKmers= new HashSet<kmerData>();
			for (kmerData kmerData : onePrefix) {
				for(seedData finalSeed : generateSeeds_extend(new seedData(kmerData), kmers)){
					usedKmers.addAll(finalSeed.getKmers());
					if(finalSeed.size()>=minSeedSize){
						seedCount++;
						out.write(">"+seedCount+" "+finalSeed.size()+" cyclic\n"+finalSeed.consensus()+"\n");
						outLong.write(">"+seedCount+" "+finalSeed.size()+" cyclic\n"+finalSeed.toString()+"\n");
					}
				}
			}
			//clean
			uniqueKmers.removeAll(usedKmers);
			
			for (kmerData kmerData : usedKmers) {
				kmers.remove(kmerData.prefix());
				kmers.remove(kmerData.prefixRevComp());
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
		
		out.close();
		outLong.close();
//		return seeds;
		System.err.println("Done!");
	}
	
	private static ArrayList<seedData> generateSeeds_extend(seedData seed, HashMap<String,HashSet<kmerData>> kmers){
//		System.err.print("\t\t");
//		for(int i=0;i<seed.size();i++){
//			System.err.print(" ");
//		}
//		System.err.print("1");
//		for(int i=seed.size()+2;i<seed.getKmers().get(0).getKmer().length();i++){
//			System.err.print(" ");
//		}
//		System.err.print("\r");
		seedData newSeed= new seedData(seed);
		ArrayList<seedData> seeds= new ArrayList<seedData>();
		boolean fin=true;
		if(kmers.containsKey(newSeed.suffix())){
			for (kmerData kmer : kmers.get(newSeed.suffix())) {
				if(newSeed.addRight(kmer)){
					seeds.addAll(generateSeeds_extend(newSeed,kmers));
					newSeed= new seedData(seed);
					fin=false;
				}
			}
		}	
//		for(int i=0;i<kmers.size();i++){
//			if(newSeed.addRight(kmers.get(i))){
//				seeds.addAll(generateSeeds_extend(newSeed,kmers));
//				newSeed= new seedData(seed);
//				fin=false;
//			}
//			if(newSeed.addLeft(kmers.get(i))){
//				seeds.addAll(generateSeeds_extend(newSeed,kmers));
//				newSeed= new seedData(seed);
//				fin=false;
//			}
//		}
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
	public static String getPrefix(String fastqFile)throws Exception{
		fastqParser fqp;
		FastqSeq fqs;
		String prefix="";
		boolean done=true;
		for (int i=6;i>-1&&done;i--){
			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),"");
			fqs=fqp.next();
			prefix=fqs.getQname().substring(0,i);
			done=true;
			for(int j=0;j<10&&fqp.hasNext();j++){
				fqs=fqp.next();
				if(!prefix.equals(fqs.getQname().subSequence(0, i))){
					done=false;
				}
			}
		}
		return prefix;
	}
}

class seedData{
	
	private ArrayList<kmerData> kmers;
	private ArrayList<Integer> strand;
	private String consensus;
	
	private seedData(){
		kmers= new ArrayList<kmerData>();
		strand= new ArrayList<Integer>();
		consensus="";
	}
	
	public seedData(kmerData first){
		this();
		kmers.add(first);
		strand.add(0);
		consensus=first.getKmer();
	}
	
	public seedData(seedData old){
		this();
		this.kmers.addAll(old.kmers);
		this.strand.addAll(old.strand);
		this.consensus=old.consensus;
	}
	
	public String suffix(int n){
//		System.err.println(n+";"+consensus.length());
		return consensus.substring(consensus.length()-n,consensus.length());
	}
	
	public String suffix(){
		return suffix(kmers.get(0).getKmer().length()-1);
	}
	
	public ArrayList<kmerData> getKmers(){
		return kmers;
	}
	
	public String consensus(){
		//presumes that the kmers are of the same length
//		String seq=kmers.get(0).getKmer();
//		int lastPos=seq.length()-1;
//		for(int i=1;i<kmers.size();i++){
//			seq+=""+kmers.get(i).getKmer().charAt(lastPos);
//		}
//		return seq;
		return consensus;
	}
	
	public boolean addRight(kmerData kd){
		if(this.contains(kd)){
			return false;
		}
		//same strand
		if(consensus.endsWith(kd.prefix())){
			kmers.add(kd);
			strand.add(0);
			consensus+=kd.tail();
			return true;
		}
		//opposite strand
		if(consensus.endsWith(kd.prefixRevComp())){
			kmers.add(kd);
			strand.add(1);
			consensus+=kmerUtils.complement(kd.head());
			return true;
		}
		return false;
	}
	
	public boolean addLeft(kmerData kd){
		if(this.contains(kd)){
			return false;
		}
		//same strand
		if(consensus.startsWith(kd.suffix())){
			kmers.add(0, kd);
			strand.add(0,0);
			consensus=kd.head()+consensus;
			return true;
		}
		//opposite strand
		if(consensus.startsWith(kd.suffixRevComp())){
			kmers.add(0, kd);
			strand.add(0,1);
			consensus=kmerUtils.complement(kd.tail())+consensus;
		}
		return false;
	}
	
//	public boolean addRight(kmerData ext){
//		if(extendRight(ext)&&!contains(ext)){
//			kmers.add(ext);
//			return true;
//		}
//		return false;
//	}
//	
//	public boolean addLeft(kmerData ext){
//		if(extendLeft(ext)&&!contains(ext)){
//			kmers.add(0, ext);
//			return true;
//		}
//		return false;
//	}
//	
//	public boolean extendRight(kmerData ext){
//		return kmers.get(kmers.size()-1).extendRight(ext);
//	}
//	
//	public boolean extendLeft(kmerData ext){
//		return kmers.get(0).extendLeft(ext);
//	}
	
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
		for (int i=0;i<kmers.size();i++){
			if(strand.get(i)==0){
				out+="\n"+kmers.get(i).toString();
			}else{
				out+="\n"+kmers.get(i).toStringRevComp();
			}
		}
		
//		String lastSuffix=kmers.get(0).suffix();
//		kmerData kmer;
//		for(int i=1;i<kmers.size();i++){
//			kmer=kmers.get(i);
//			if(kmer.getKmer().startsWith(lastSuffix)){
//				out+="\n"+kmer.toString();
//			}else{
//				out+="\n"+kmer.toStringRevComp();
//			}
//			lastSuffix=kmer.suffix();
//		}
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
	private String kmerSeqRevComp;
	private String[] data;
	
	public kmerData(int pos, String s) throws Exception{
		this(pos, s.split("\t"));
	}
	
	public kmerData(int pos, String[] data) throws Exception {
		super();
		this.pos = pos;
		this.data = data;
		kmerSeqRevComp=kmerUtils.reverseComplementSeq(data[pos]);
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
	
	public kmerData revComp()throws Exception{
		return new kmerData(pos, this.toStringRevComp());
	}
	
	public kmerData reduce(int n) throws Exception{
		String[] newData=new String[data.length];
		for(int i=0;i<data.length;i++){
			newData[i]=data[i];
		}
		newData[pos]=kmerUtils.kmerToUse(data[pos]);
		return new kmerData(pos,newData);
	}
	
	public String suffixRevComp(){
		return suffixRevComp(1);
	}
	
	public String suffixRevComp(int n){
		return getKmerRevComp().substring(n);
	}
	
	public String prefixRevComp(int n){
		return getKmerRevComp().substring(0, getKmer().length()-n);
	}
	
	public String prefixRevComp(){
		return prefixRevComp(1);
	}
	
	public String tail(){
		return ""+getKmer().charAt(getKmer().length()-1);
	}
	
	public String head(){
		return ""+getKmer().charAt(0);
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
	
	private boolean extendLeft(String prefix){
		return getKmer().startsWith(prefix)||getKmerRevComp().startsWith(prefix);
	}
	
	public boolean extendLeft(kmerData prefix,int n){
		return extendLeft(prefix.suffix(n));
//		return startsWith(prefix.getKmer().substring(n));
	}
	
	public boolean extendLeft(kmerData prefix){
		return extendLeft(prefix,1);
	}
	
	private boolean extendRight(String suffix){
		return getKmer().endsWith(suffix)||getKmerRevComp().endsWith(suffix);
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
	
	private String getKmerRevComp(){
		return kmerSeqRevComp;
	}

	public int compareTo(kmerData o) {
		if(this.rep.equals(o.rep))
			return this.permutations-o.permutations;
		return this.rep.compareTo(o.rep);
	}
	
	public String toString(){
		return this.toString("\t");
	}
	
	public String toStringRevComp(){
		return toStringRevComp("\t");
	}
	
	public String toStringRevComp(String sep){
		String s="";
		if(pos==0){
			s=getKmerRevComp();
		}else{
			s=data[0];
		}
		for(int i=1;i<data.length;i++){
			if(i==pos){
				s+=sep+getKmerRevComp();
			}else{
				s+=sep+data[i];
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
