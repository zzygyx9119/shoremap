package tools.fastq;

import information.StringContent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.regex.MatchResult;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import tools.fasta.FastaSeq;
import tools.fasta.fastaParser;
import tools.kmer.KmerSet_binary;
import tools.kmer.KmerTable;
import tools.kmer.Kmer_count_tuple;
import tools.reptile.ReptileLine;
import tools.rocheQual.RocheQualParser;
import tools.rocheQual.RocheQualSeq;

public class fastqUtils {

	private static final String sep="\t";
	
	public static void main(String[] args) throws Exception{
		if(args.length>0){
			if(args[0].equals("information")&&args.length>7){
				information(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]),Integer.parseInt(args[4]),Integer.parseInt(args[5]),Integer.parseInt(args[6]),args[7]);
			}else if(args[0].equals("tabulatedTopOnly")&&args.length==2){
//				tabulatedTopOnly(args[1]);
			}else if(args[0].equals("sub")&&args.length>2){
				ArrayList<String> fastqFiles= new ArrayList<String>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(args[i]);
				}
				sub(fastqFiles, args[args.length-1], true);
			}else if(args[0].equals("subJia")&&args.length>3){
				subJia(args[1],args[2],true,Integer.parseInt(args[3]));
			}else if(args[0].equals("antisub")&&args.length>2){
				ArrayList<String> fastqFiles= new ArrayList<String>();
				for(int i=1;i<args.length-1;i++){
					fastqFiles.add(args[i]);
				}
				sub(fastqFiles, args[args.length-1], false);
			}else if(args[0].equals("subReg")&&args.length>2){
				subReg(args[1],args[2],true);
			}else if(args[0].equals("antisubReg")&&args.length>2){
				subReg(args[1],args[2],false);
			}else if(args[0].equals("clean")&&args.length>6){
				clean(args[1],args[2],Integer.parseInt(args[3]),Integer.parseInt(args[4]),args[5],Integer.parseInt(args[6]));
			}else if(args[0].equals("randomsub")&&args.length>2){
				randomsub(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("reverseComplement")&&args.length>1){
				reverseComplement(args[1]);
			}else if(args[0].equals("removeShort")&&args.length>2){
				removeShort(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("removeShortVelvetPair")&&args.length>4){
				removeShortVelvetPair(args[1],Integer.parseInt(args[2]),args[3],args[4]);
			}else if(args[0].equals("convertNumberToPhredQual")&&args.length>1){
				convertNumberToPhredQual(args[1]);
			}else if(args[0].equals("randomVelvetPairedSubset")&&args.length>2){
				randomVelvetPairedSubset(args[1],Double.parseDouble(args[2]));
			}else if(args[0].equals("kmerCount")&&args.length>4){
				kmerCount(args[1],Integer.parseInt(args[2]),args[3],args[4]);
			}else if(args[0].equals("kmerCount_histogram")&&args.length>2){
				kmerCount_histogram(args[1], Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerCount_merge")&&args.length>3){
				ArrayList<String> files= new ArrayList<String>();
				for(int i=1;i<args.length-1;i++){
					files.add(args[i]);
				}
				kmerCount_merge(files,args[args.length-1]);
			}else if(args[0].equals("kmerize")&&args.length>2){
				kmerizeEnd(args[1],Integer.parseInt(args[2]),Integer.parseInt(args[3]));
			}else if(args[0].equals("kmerizeEnd")&&args.length>3){
				kmerize(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerCoverage")&&args.length>3){
				ArrayList<String> fastqFiles=new ArrayList<String>();
				for(int i=3;i<args.length;i++){
					fastqFiles.add(args[i]);
				}
				kmerCoverage(fastqFiles, readKmers(args[1], Integer.parseInt(args[2])));
			}else if(args[0].equals("kmerCoverage_count")&&args.length>3){
				ArrayList<String> fastqFiles=new ArrayList<String>();
				for(int i=3;i<args.length;i++){
					fastqFiles.add(args[i]);
				}
				kmerCoverage_count(fastqFiles, readKmers(args[1], Integer.parseInt(args[2])));
			}else if(args[0].equals("kmerCoverage_stored")&&args.length>2){
				ArrayList<String> fastqFiles=new ArrayList<String>();
				for(int i=2;i<args.length;i++){
					fastqFiles.add(args[i]);
				}
				kmerCoverage(fastqFiles, readKmers(args[1]));
			}else if(args[0].equals("storeKmerSet")&&args.length>3){
				storeKmerSet(args[1],Integer.parseInt(args[2]),args[3]);
			}else if(args[0].equals("kmerCoverage_generateCorrection")&&args.length>2){
				kmerCoverage_generateCorrection(args[1],Integer.parseInt(args[2]));
			}else if(args[0].equals("kmerCoverage_extractPerfect")&&args.length>3){
				kmerCoverage_extractPerfect(args[1],args[2],Integer.parseInt(args[3]));
			}else if(args[0].equals("reptileCorrection")&&args.length>2){
				reptileCorrection(args[1],args[2]);
			}else if(args[0].equals("benchmark")&&args.length==3){
				benchmark(args[1],args[2]);
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
		String help="Usage: illuminaUtils <cmd> <input>\n";
		help+="where <cmd> is:\n";
		help+="information - \n";
		help+="\t<input> = <fastqFile> <burnin> <max nr to sample> <sample  distance> <min word length (inclusive)> <max word length (exclusive)> <outPrefix>\n";
		help+="randomsub - extracts a random subset from a fastq file. Each sequence will be given one chance in N.\n";
		help+="\t<input> = <fastqFile> <N>\n";
		help+="(anti)sub - extracts a subset from a fastq file\n";
		help+="\t<input> = <fastqFile> <subset file>\n";
		help+="(anti)subReg - extracts a subset from a fastq file\n";
		help+="\t<input> = <fastqFile> <regular expression>\n";
		help+="clean - cleans a fastq file from stretches of the letters in the word. That is, XAT will search for stretches X's, A's or T's with a minimum length of word size in the fasta file and return the reads from the fastq file. Also ditches reads with too many N's\n";
		help+="\t<input> = <fastqFile> <fastaFile> <min length> <word size> <word> <max number of N's>\n";
		help+="reverseComplement - reverse complements a fastqfile\n";
		help+="\t<input> = <fastqFile> \n";
		help+="removeShort - drop sequences shorter than n\n";
		help+="\t<input> = <fastqFile> <n>\n";
		help+="removeShortVelvetPair - drop sequences shorter than n\n";
		help+="\t<input> = <fastqFile> <n> <singletFile> <pairFile>\n";
		help+="convertNumberToPhredQual - converts number qualities to Phred values\n";
		help+="\t<input> = <fastqFile>\n";
		help+="randomVelvetPairedSubset - Selects a random number of reads corresponding to the given fraction (0-1.0) from a paired (interleaved velvet style) fastq file \n";
		help+="\t<input> = <fastqFile> <fraction>\n";
//		help+="kmerCount - counts the kmers in a fastqfile and writes a java object file\n";
//		help+="\t<input> = <fastqFile> <kmerSize> <tmpPrefix> <outFile>\n";
//		help+="kmerCount_merge - merges several kmerize|sort|uniq -c|awk '{print $1\"\\t\"$2}' files\n";
//		help+="\t<input> = <kmer uniq file 1> <...> <kmer uniq file n> <outFile>\n";
		
		
		
//		help+="kmerCount_histogram - makes a histogram over the kmerCounts\n";
//		help+="\t<input> = <kmer obj file> <bin size>\n";
		
		help+="kmerize - prints all kmers of size n in each sequence\n";
		help+="\t<input> = <fastqFile> <n>\n";
		help+="kmerizeEnd - prints the kmers in the end of a read missing when reducing from n to m\n";
		help+="\t<input> = <fastqFile> <n> <m>\n";
		help+="kmerCoverage - prints kmer-coverage of all sequences in a file with respect to a KmerSet file\n";
		help+="\t<input> = <KmerCountFile> <cutoff> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerCoverage_stored - prints kmer-coverage of all sequences in a file with respect to a KmerSet file\n";
		help+="\t<input> = <KmerCountFile> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerCoverage_count - prints total kmer-coverage (count) of all sequences in a file with respect to a KmerSet file\n";
		help+="\t<input> = <KmerCountFile> <cutoff> <fastqFile 1> ... <fastqFile n>\n";
		help+="kmerCoverage_generateCorrection - reads the output from kmerCoverage and generates a file with the longest fully covered subsequence\n";
		help+="\t<input> = <KmerCoverage file> <kmerSize>\n";
		help+="kmerCoverage_extractPerfect - prints the qname of all sequences with perfect kmer-coverage\n";
		help+="\t<input> = <fastqFile> <KmerCoverage file> <kmerSize>\n";
		help+="storeKmerSet - stores the kmers (count\\tkmer) that occcurs at least <min cutoff> times in a kmerSet on the disk for future use\n";
		help+="\t<input> = <Kmer count file> <min cutoff> <outfile>\n";
		help+="reptileCorrection - Applies the corrections calculated by reptile (Assumes that the read index starts at 0)\n";
		help+="\t<input> = <fastqFile> <reptileFile>\n";
		
		
		
		return help;
	}
	
	public static void benchmark(String fastqNew,String fastqOld)throws Exception{
		FastqSeq fqs;
		long start= System.currentTimeMillis();
		FastqParser_nio fqpN1= new FastqParser_nio(fastqNew, getPrefix(fastqNew));
		for(;fqpN1.hasNext();){
			fqs=fqpN1.next();
		}
		System.out.println("First file, nio-parser: "+(System.currentTimeMillis()-start));
		fqpN1.printTimes();
		start=System.currentTimeMillis();
		fastqParser fqpO1= new fastqParser(new BufferedReader(new FileReader(fastqNew)),getPrefix(fastqNew));
		for(;fqpO1.hasNext();){
			fqs=fqpO1.next();
		}
		System.out.println("First file, old parser: "+(System.currentTimeMillis()-start));
		start=System.currentTimeMillis();
		fastqParser fqpO2= new fastqParser(new BufferedReader(new FileReader(fastqOld)), getPrefix(fastqOld));
		for(;fqpO2.hasNext();){
			fqs=fqpO2.next();
		}
		System.out.println("Second file, old parser: "+(System.currentTimeMillis()-start));
		start=System.currentTimeMillis();
		FastqParser_nio fqpN2= new FastqParser_nio(fastqOld, getPrefix(fastqOld));
		for(;fqpN2.hasNext();){
			fqs=fqpN2.next();
		}
		System.out.println("Second file, nio-parser: "+(System.currentTimeMillis()-start));
		fqpN2.printTimes();
	}
	
	public static void kmerCoverage_extractPerfect(String fastqFile, String coverageFile,int kmerSize)throws Exception{
		RocheQualParser rqp= new RocheQualParser(new BufferedReader(new FileReader(coverageFile)));
		RocheQualSeq rqs;
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		String[] quals;
		boolean perfect;
		for(int i=0;rqp.hasNext()&&fqp.hasNext();i++){
			if(i%10000==0){
				System.err.print("\t"+i+"\r");
			}
			rqs=rqp.next();
			
			quals=rqs.getSeq().trim().split(" +");
			perfect=true;
			if(Integer.parseInt(quals[0])==1&&Integer.parseInt(quals[quals.length-1])==1&&Integer.parseInt(quals[quals.length-kmerSize])==kmerSize){
				for(int j=kmerSize;j<quals.length-kmerSize;j+=kmerSize-1){
					if(Integer.parseInt(quals[j])!=kmerSize){
						perfect=false;
						break;
					}
				}
				if(perfect){
					System.out.println(fqp.next());
				}else{
					fqp.next();
				}
			}else{
				fqs=fqp.next();
				if(!rqs.getQname().equals(fqs.getQname())){
					System.err.println("Out of sync, q: "+rqs.getQname()+", fq: "+fqs.getQname());
				}
			}

		}
	}

	public static void kmerCoverage_generateCorrection(String coverageFile,int kmerSize)throws Exception{
		RocheQualParser rqp= new RocheQualParser(new BufferedReader(new FileReader(coverageFile)));
		RocheQualSeq rqs;
		int[] coords,quals;
		for(;rqp.hasNext();){
			rqs=rqp.next();
			quals=rqs.getIntSeq();
			coords=kmerCoverage_generateCorrection_longestMaxSub(quals);
			System.out.println(rqs.getQname()+sep+quals[coords[0]]+sep+coords[0]+sep+coords[1]+sep+kmerCoverage_count(quals,kmerSize));
		}
	}
	
	public static int kmerCoverage_count(int[] quals,int kmerSize){
		int count=quals[0];
		for(int i=1;i<quals.length;i++){
			if(quals[i]>quals[i-1]||quals[i]==kmerSize){
				count++;
			}
		}
		return count;
	}
	
	
	public static int[] kmerCoverage_generateCorrection_longestMaxSub(int[] quals){
		int[] coords= new int[]{0,0};
		int[] tmpcoords=new int[]{-1,-1};
		for(int i=0;i<quals.length;i++){
			if(quals[i]>quals[coords[0]]){
				//found a new max
				coords[0]=i;
				coords[1]=i+1;
				tmpcoords=new int[]{-1,-1};
			}else if(quals[i]==quals[coords[0]]){
				if(i==coords[1]){
					coords[1]++;
				}else if(tmpcoords[0]!=-1){
					//start a second stretch
					tmpcoords[0]=i;
					tmpcoords[1]=i+1;
				}else if(tmpcoords[1]==i){
					tmpcoords[1]++;
				}else{
					//check if a longer stretch is found and start a new
					if(tmpcoords[1]-tmpcoords[0]>coords[1]-coords[0]){
						coords[0]=tmpcoords[0];
						coords[1]=tmpcoords[1];
					}
					tmpcoords[0]=i;
					tmpcoords[1]=i+1;
				}
			}
		}
		if(tmpcoords[1]-tmpcoords[0]>coords[1]-coords[0]){
			coords[0]=tmpcoords[0];
			coords[1]=tmpcoords[1];
		}
		return coords;
	}
	
	public static void reptileCorrection(String fastqFile,String reptileFile)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		int nextToBeRead=0,changed=0,base=0;
		ReptileLine rl;
		BufferedReader reptileIn= new BufferedReader(new FileReader(reptileFile));
		for(String s=reptileIn.readLine();s!=null;s=reptileIn.readLine()){
			rl=new ReptileLine(s);
			for(;fqp.hasNext()&&nextToBeRead<rl.getReadId()-base;){
				System.out.println(fqp.next());
				nextToBeRead++;
			}
			if(fqp.hasNext()){
				fqs=fqp.next();
				nextToBeRead++;
				fqs.correct(rl.getCorrections());
				changed++;
				System.out.println(fqs);
			}else{
				throw new Exception("No reads to cover for error with readID: "+rl.getReadId());
			}
		}
		for(;fqp.hasNext();){
			System.out.println(fqp.next());
		}
		System.err.println("Corrected "+changed+" sequences");
	}
	
	private static KmerSet_binary readKmers(String kmerFile)throws Exception{
		return (KmerSet_binary)(new ObjectInputStream(new FileInputStream(kmerFile))).readObject();
	}
	
	public static void kmerCoverage_count(ArrayList<String> fastqFiles,KmerSet_binary kmers) throws Exception{
		fastqParser fqp;
		FastqSeq fqs;
		BufferedWriter out;
		for (String fastqFile : fastqFiles) {
			System.err.println(fastqFile);
			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
			out= new BufferedWriter(new FileWriter(fastqFile+".kmerCoverageCount"));
			for(int i=0;fqp.hasNext();i++){
				if(i%10000==0){
					System.err.print("\t"+i+"\r");
				}
				fqs=fqp.next();
				out.write(fqs.getQname()+"\t"+kmers.count(fqs.getSeq())+"\n");
			}
		}
	}
	
	public static void kmerCoverage(ArrayList<String> fastqFiles,KmerSet_binary kmers)throws Exception{
		System.err.println("Reading background file...");
//		KmerSet kmers= (KmerSet)(new ObjectInputStream(new FileInputStream(kmerBackgroundFile))).readObject();
		System.err.println("Parsing sequence file...");
		fastqParser fqp;
		FastqSeq fqs;
		int[] cov;
		BufferedWriter out;
		for (String fastqFile : fastqFiles) {
			System.err.println(fastqFile);
			fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
			out= new BufferedWriter(new FileWriter(fastqFile+".kmerCoverage"));
			for(int i=0;fqp.hasNext();i++){
				if(i%10000==0){
					System.err.print("\t"+i+"\r");
				}
				fqs=fqp.next();
				out.write(">"+fqs.getQname()+"\n");
				cov=kmers.coverage(fqs.getSeq());
				out.write(""+cov[0]);
				for(int j=1;j<cov.length;j++){
					out.write(" "+cov[j]);
				}
				out.write("\n");
			}
			out.close();
		}
		System.err.println("\nDone!");
	}
	
	public static int[] kmerCoverage_coverage(String seq,KmerSet_binary kmers){
		int kmerSize=kmers.getKmerSize();
		int[] cov=new int[seq.length()];
		boolean[] goodKmer= new boolean[seq.length()-kmerSize+1];
		String kmer=seq.substring(0, kmerSize);
		goodKmer[0]=kmers.exists(kmer);
		cov[0]=goodKmer[0]?1:0;
		for(int i=1,j=kmerSize;j<seq.length();i++,j++){
			kmer=kmer.substring(1)+seq.charAt(j);
			goodKmer[i]=kmers.exists(kmer);
			cov[i]=cov[i-1]+(goodKmer[i]?1:0)-(i>=kmerSize?(goodKmer[i-kmerSize]?1:0):0);
		}
		for(int i=seq.length()-kmerSize+1;i<seq.length();i++){
			cov[i]=cov[i-1]-(goodKmer[i-kmerSize]?1:0);
		}
		return cov;
	}
	
	private static KmerSet_binary readKmers(String kmerCountFile, int cutoff)throws Exception{
		BufferedReader in= new BufferedReader(new FileReader(kmerCountFile));
		String[] l=in.readLine().split("\t");
		KmerSet_binary kmers= new KmerSet_binary(l[1].length());
		if(Integer.parseInt(l[0])>=cutoff){
			kmers.addSeq(l[1]);
		}
		int line=1;
		for(String s=in.readLine();s!=null;s=in.readLine(),line++){
			if(line%10000==0)
				System.err.print("\t"+line+"   "+kmers.size()+"\r");
			l=s.split("\t");
			if(Integer.parseInt(l[0])>=cutoff){
				kmers.addSeq(l[1]);
			}
		}
		
		return kmers;
	}
	
	public static void storeKmerSet(String kmerCountFile,int cutoff,String outfile)throws Exception{
		KmerSet_binary kmers= readKmers(kmerCountFile, cutoff);
		System.err.println("\nwriting to file...");
		ObjectOutputStream out= new ObjectOutputStream(new FileOutputStream(outfile));
		out.writeObject(kmers);
		out.close();
	}
	
	public static void kmerize(String fastqFile, int n)throws Exception{
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		for(;fqp.hasNext();){
			String seq=fqp.next().getSeq();
			for(int i=0,j=n;j<=seq.length();i++,j++){
				System.out.println(seq.substring(i, j));
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
				System.out.println(seq.substring(i, j));
			}
		}
	}
	
	public static void kmerCount(String fastqFile, int kmerSize,String tmpPrefix,String outFile)throws Exception{
		fastqParser fqp=new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		KmerTable kTable= new KmerTable(kmerSize,tmpPrefix,1000000);
		for(int i=1;fqp.hasNext();i++){
			if(i%10000==0){
//				System.err.print("\t"+i+"\t"+kTable.size()+"\r");
			}
			kTable.addSequence(fqp.next().getSeq());
		}
		System.err.println();
		kTable.merge(outFile);
//		kTable.printHistogram();
	}
	
	public static void kmerCount_merge(ArrayList<String> files,String outFile)throws Exception{
		KmerTable kTable= new KmerTable(files);
		kTable.merge(outFile);
	}
	
	public static void kmerCount_histogram(String kmerFile,int bin)throws Exception{
		ObjectInputStream in= new ObjectInputStream(new FileInputStream(kmerFile));
		HashMap<Integer, Integer>counts= new HashMap<Integer, Integer>();
		int max=0,curBin;
		Kmer_count_tuple kmer;
		try{
			for(int i=0;true;i++){
				if(i%10000==0){
					System.err.print("\t"+i+"\r");
				}
				kmer=(Kmer_count_tuple) in.readObject();
				curBin=kmer.getCount()/bin;
				max=curBin>max?curBin:max;
				if(counts.containsKey(curBin)){
					counts.put(curBin, counts.get(curBin)+1);
				}else{
					counts.put(curBin, 1);
				}
			}
		}catch (EOFException e) {
		}
		System.err.println("printing...");
		//print
		System.out.println("kmer count (bin size="+bin+")\t#unique kmers");
		for(int i=1;i<=max;i++){
			if(counts.containsKey(i)){
				System.out.println(i*bin+"\t"+counts.get(i));
			}else{
				System.out.println(i*bin+"\t0");
			}
		}
	}
	
	public static void randomVelvetPairedSubset(String fastqFile,double fraction)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		Random rg= new Random();
		for(;fqp.hasNext();){
			if(rg.nextDouble()<fraction){
				System.out.println(fqp.next());
				System.out.println(fqp.next());
			}else{
				fqp.next();
				fqp.next();
			}
		}
	}
	
	private static void convertNumberToPhredQual(String fastqFile)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		for(;fqp.hasNext();){
			fqs=fqp.next();
			fqs.convertNumberToPhredQual();
			System.out.println(fqs.toString());
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
	
	private static void removeShortVelvetPair(String fastqFile, int n,String singletFile,String pairFile)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs1,fqs2;
		BufferedWriter singlet= new BufferedWriter(new FileWriter(singletFile));
		BufferedWriter pair= new BufferedWriter(new FileWriter(pairFile));
		for(;fqp.hasNext();){
			fqs1=fqp.next();
			if(fqp.hasNext()){
				fqs2=fqp.next();
				if(fqs1.length()<n){
					if(fqs2.length()>=n){
						singlet.write(fqs2.toString()+"\n");
					}
				}else{
					if(fqs2.length()<n){
						singlet.write(fqs1.toString()+"\n");
					}else{
						pair.write(fqs1.toString()+"\n");
						pair.write(fqs2.toString()+"\n");
					}
				}
			}else{
				throw new Exception("File does not contain an even number of reads... last read: "+fqs1.getQname());
			}
		}
		singlet.close();
		pair.close();
	}
	
	private static void removeShort(String fastqFile, int n)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		for(;fqp.hasNext();){
			fqs= fqp.next();
			if(fqs.length()>=n){
				System.out.println(fqs.toString());
			}
		}
	}
	
	private static void reverseComplement(String fastqFile)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		for(;fqp.hasNext();){
			System.out.println(fqp.next().reverseComplement().toString());
		}
	}
	
	private static void clean(String fastqFile,String fastaFile,int minLength,int wordSize, String word, int maxNs) throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		fastaParser fp= new fastaParser(new BufferedReader(new FileReader(fastaFile)));
		FastqSeq fqs;
		FastaSeq fs;
//		ArrayList<Pattern> toMask= new ArrayList<Pattern>();
		ArrayList<Matcher> Masks= new ArrayList<Matcher>();
		
		for(int i=0;i<word.length();i++){
			Masks.add(Pattern.compile(word.charAt(i)+"{"+wordSize+",}").matcher(""));
		}
		MatchResult mr;
		ArrayList<int[]> toMask;
		Comparator<int[]> pairSorter= new IntPairComparator();
		int printed,start,end,nCount;
		//assumes all reads are present in both files
		for(;fqp.hasNext()&&fp.hasNext();){
			fqs=fqp.next();
			fs=fp.next();
			if(!fs.getQname().equals(fqs.getQname())){
				System.err.println("Unmatched names - fasta: "+fs.getQname()+", fastq: "+fqs.getQname());
			}
			nCount=0;
			for (int i=0;i<fs.getSeq().length();i++) {
				if(fs.getSeq().charAt(i)=='N'){
					nCount++;
				}
			}
			if(nCount<=maxNs){
				toMask= new ArrayList<int[]>();
				toMask.add(new int[]{0,0});
				toMask.add(new int[]{fs.length(),fs.length()});
				for (Matcher matcher : Masks) {
					matcher.reset(fs.getSeq());
					while(matcher.find()){
						mr=matcher.toMatchResult();
						toMask.add(new int[]{mr.start(),mr.end()});
					}
				}
				Collections.sort(toMask, pairSorter);
				printed=0;
				for(int i=1;i<toMask.size();i++){
					start=toMask.get(i-1)[1];
					end=toMask.get(i)[0];
					if(end-start>minLength){
						//print
						printed++;
						if(printed==1){
							//without change to name
							System.out.println(fqs.toString(start, end));
						}else{
							//with added idenifier
							System.out.println(fqs.toString(start, end, fqs.getQname()+"_"+printed));
						}
					}
				}
			}
		}
	}
	
	private static void randomsub(String fastqFile,int N)throws Exception{
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqParser_nio fqp= new FastqParser_nio(fastqFile, getPrefix(fastqFile));
		Random rs= new Random();
		final int luckynumber=rs.nextInt(N);
		for (;fqp.hasNext();){
			if(rs.nextInt(N)==luckynumber){
				System.out.println(fqp.next().toString());
			}else{
				fqp.next();
			}
		}
	}
	
	private static void sub(ArrayList<String> fastqFiles,String subsetFile,boolean inSubset) throws Exception{
		HashSet<String> subset= fileToHashset(subsetFile);
		for (String fastqFile : fastqFiles) {
			sub(fastqFile,subset,inSubset);	
		}
	}
	
//	private static void sub(String fastqFile,String subsetFile, boolean inSubset)throws Exception{
//		sub(fastqFile,fileToHashset(subsetFile),inSubset);
		
//		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
//		FastqSeq fqs;
//		
//		for(;fqp.hasNext();){
//			fqs=fqp.next();
//			if(subset.contains(fqs.getQname())==inSubset){
//				System.out.println(fqs.toString());
//			}
//		}
//	}
	
	private static HashSet<String> fileToHashset(String file)throws Exception{
		HashSet<String> hashset= new HashSet<String>();
		BufferedReader in= new BufferedReader(new FileReader(file));
		for(String s= in.readLine();s!=null;s=in.readLine()){
			hashset.add(s);
		}
		return hashset;
	}
	
	private static void sub(String fastqFile, HashSet<String> subset, boolean inSubset)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		
		for(;fqp.hasNext();){
			fqs=fqp.next();
			if(subset.contains(fqs.getQname())==inSubset){
				System.out.println(fqs.toString());
			}
		}
	}
	
	private static void subJia(String fastqFile,String subsetFile, boolean inSubset,int dropPost)throws Exception{
		HashSet<String> subset= fileToHashset(subsetFile);
		
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		
		for(;fqp.hasNext();){
			fqs=fqp.next();
			if(subset.contains(fqs.getQname().substring(0, fqs.getQname().length()-dropPost))==inSubset){
				System.out.println(fqs.toString());
			}
		}
	}
	
	private static void subReg(String fastqFile, String regExpString, boolean inSubset) throws Exception{
		Matcher regExpMatcher=Pattern.compile(regExpString).matcher("");
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		for(;fqp.hasNext();){
			fqs=fqp.next();
			regExpMatcher.reset(fqs.getQname());
			if(regExpMatcher.find()==inSubset){
				System.out.println(fqs.toString());
			}
		}
	}
	
	private static void information(String fastqFile,int burnin,int sampleSize,int sampleDistance,int minWordLength, int maxWordLength,String outPrefix)throws Exception{
		fastqParser fqp= new fastqParser(new BufferedReader(new FileReader(fastqFile)),getPrefix(fastqFile));
		FastqSeq fqs;
		StringContent seqInfo= new StringContent();
//		StringContent qual= new StringContent(minWordLength);
		final int lastSample= sampleSize==0?Integer.MAX_VALUE:burnin+sampleSize*sampleDistance+1; 
		for(int count=1;fqp.hasNext()&&count<lastSample;count++){
			if(count>burnin && (count-burnin)%sampleDistance==0){
				//sample
				fqs=fqp.next();
				final String seq=fqs.getSeq();
				final int length=seq.length();
				for(int wordLength=minWordLength;wordLength<maxWordLength;wordLength++){
					for(int i=0;i<length-wordLength;i++){
						seqInfo.add(seq.substring(i, i+wordLength));
					}
				}
			}else{
				//discard
				fqp.next();
			}
		}
		//print
		for(int wordLength=minWordLength;wordLength<maxWordLength;wordLength++){
			BufferedWriter out= new BufferedWriter(new FileWriter(outPrefix+"_wordLength_"+wordLength+".csv"));
			out.write("#Total"+sep+seqInfo.getTotal(wordLength+"")+"\nWord"+sep+"Count\n");
			HashMap<String, Integer> dist=seqInfo.getCount(wordLength+"");
			for(String key : dist.keySet()){
				out.write(key+sep+dist.get(key)+"\n");
			}
			out.close();
		}
	}
}

class IntPairComparator implements Comparator<int[]>{

	public int compare(int[] arg0, int[] arg1) {
		//assumes that the arrays will be of size two, but at least one.
		return arg0[0]-arg1[0];
	}
	
}
