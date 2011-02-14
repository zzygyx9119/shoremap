package tools.kmer;

import java.util.HashMap;

import utils.BloomFilter;

public class KmerTable_hash {

	private HashMap<boolean[], Integer> counts;
	private BloomFilter<boolean[]> preCheck;
	private int kmerSize,prefixSize;
	private String tmpPrefix;
	private KmerTable_utils kt_utils;
	
	
	public KmerTable_hash(){
		this(31);
	}
	
	public KmerTable_hash(int kmerSize){
		this(kmerSize,"tmp_kmerTable",4);
	}
	
	public KmerTable_hash(int kmerSize, String tmpPrefix,int prefixSize){
		
//		if(new Kmer("TTAGACCGAATGATCTTGTTTTCATCCCTTC").equals(new Kmer("TTAGACCGAATGATCTTGTTTTCATCCCTTC"))){
//			System.out.println("tjo");
//		}
		this.prefixSize=prefixSize;
		this.kmerSize=kmerSize;
		this.tmpPrefix=tmpPrefix;
		kt_utils=new KmerTable_utils(kmerSize);
		System.err.println("initializing hash...");
		counts=new HashMap<boolean[], Integer>(1000000000);
		System.err.println("initializing bloom filter...");
		preCheck=new BloomFilter<boolean[]>(10000000, (int)Math.min(Math.pow(5.0, kmerSize), Integer.MAX_VALUE));
		System.err.println("initialization done...");
	}
	
	public int getCount(String kmer_seq){
		boolean[] kmer= kt_utils.stringToBoolean(kmer_seq);
		if(exist(kmer)){
			return counts.get(kmer);
		}else{
			return 0;
		}
	}
	
	public int addSequence(String seq)throws Exception{
		if(seq.length()>kmerSize){
			boolean[] kmer= addKmer(seq.substring(0, kmerSize));
			addKmer(kmer);
			for(int i=kmerSize;i<seq.length();i++){
//				System.out.println(kmer_seq);
				kmer=addKmer(kt_utils.shiftBoolean(kmer, seq.charAt(i)));
			}
		}
		return seq.length()-kmerSize+1;
	}
	

	
	private boolean[] addKmer(String kmer_seq)throws Exception{
		//change string to something with a more optimal hash function
		boolean[] kmer=kt_utils.stringToBoolean(kmer_seq);
		return addKmer(kmer,1);
	}
	
	private boolean[] addKmer(boolean[] kmer)throws Exception{
		return addKmer(kmer, 1);
	}
	
	private boolean[] addKmer(boolean[] kmer,int count)throws Exception{
		if(exist(kmer)){
			counts.put(kmer, counts.get(kmer)+count);
		}else{
			preCheck.add(kmer);
			counts.put(kmer, count);
		}
		return kmer;
	}
	
	private boolean exist(boolean[] kmer){
		//implement bloom filter
		return counts.containsKey(kmer);
	}
	
	public void printHistogram(){
		printHistogram(1);
	}
	
	public void printHistogram(int binSize){
		HashMap<Integer, Integer> histogram= new HashMap<Integer, Integer>();
		int cur,max=Integer.MIN_VALUE;
		//collect data
		for (Integer count : counts.values()) {
//			System.err.print(count);
			cur=count/binSize;
//			System.err.println("\t"+cur);
			max=cur>max?cur:max;
			if(histogram.containsKey(cur)){
				histogram.put(cur, histogram.get(cur)+1);
			}else{
				histogram.put(cur, 1);
			}
		}
		//print
		for(int i=1;i<=max;i++){
			if(histogram.containsKey(i)){
				System.out.println(i*binSize+"\t"+histogram.get(i));
			}else{
				System.out.println(i*binSize+"\t0");
			}
		}
	}
	
	public int size(){
		return counts.size();
	}
}
