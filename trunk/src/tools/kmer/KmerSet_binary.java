package tools.kmer;

import java.io.Serializable;
//import java.io.UnsupportedEncodingException;
import java.util.BitSet;
import java.util.HashSet;

//import utils.BloomFilter;

public class KmerSet_binary implements Serializable{


	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
//	private BloomFilter<BitSet> preFilter;
	private HashSet<BitSet> kmers;
	private int kmerSize;
	private KmerSet_binary_utils kt_utils;
	
//	public KmerSet_binary(){
//		this(-1);
//	}
	
	public KmerSet_binary(int kmerSize){
		this.kmerSize=kmerSize;
		kt_utils=new KmerSet_binary_utils(kmerSize);
//		preFilter= new BloomFilter<BitSet>(500000000, 1000000000);
		kmers=new HashSet<BitSet>(1000000000);
	}
	
	public int addSeq(String seq){
		int added=0;
		BitSet kmer=kt_utils.stringToBitSet(seq.substring(0, kmerSize));
		if(addKmer(kmer)){
			added++;
		}
		for(int i=kmerSize;i<seq.length();i++){
			kmer=kt_utils.shiftBitSet(kmer, seq.charAt(i));
			if(addKmer(kmer)){
				added++;
			}
		}
		return added;
	}
	
	public boolean addKmer(String kmer_seq){
		if(kmer_seq.length()==kmerSize){
			return this.addKmer(kt_utils.stringToBitSet(kmer_seq));
		}else{
			System.err.println("Tried to add a kmer of a different length than expected:");
			System.err.println(kmer_seq);
			return false;
		}
	}
	
	private boolean addKmer(BitSet kmer){
//		try{
//			preFilter.add(kmer);
//		}catch (UnsupportedEncodingException e) {
//
//		}
		return kmers.add(kmer);
	}
	
	public int[] coverage(String seq){
		int[] cov=new int[seq.length()];
		boolean[] goodKmer= new boolean[seq.length()-kmerSize+1];
		BitSet kmer=kt_utils.stringToBitSet(seq.substring(0, kmerSize));
		goodKmer[0]=exists(kmer);
		cov[0]=goodKmer[0]?1:0;
		for(int i=1,j=kmerSize;j<seq.length();i++,j++){
			kmer=kt_utils.shiftBitSet(kmer, seq.charAt(j));
			goodKmer[i]=exists(kmer);
			cov[i]=cov[i-1]+(goodKmer[i]?1:0)-(i>=kmerSize?(goodKmer[i-kmerSize]?1:0):0);
		}
		for(int i=seq.length()-kmerSize+1;i<seq.length();i++){
			cov[i]=cov[i-1]-(goodKmer[i-kmerSize]?1:0);
		}
		return cov;
	}
	
	public int count(String seq){
		BitSet kmer=kt_utils.stringToBitSet(seq.substring(0, kmerSize));
		int count=exists(kmer)?1:0;
		for(int i=1,j=kmerSize;j<seq.length();i++,j++){
			kmer=kt_utils.shiftBitSet(kmer, seq.charAt(j));
			if(exists(kmer)){
				count++;
			}
		}
		return count;
	}
	
	private boolean exists(BitSet kmer){
//		System.out.println(kmer);
		return kmers.contains(kmer);
	}
	
	public boolean exists(String kmer_seq){
		if(kmer_seq.length()!=kmerSize){
			return false;
		}
		BitSet kmer= kt_utils.stringToBitSet(kmer_seq);
//		try{
//			if(preFilter.contains(kmer)){
//				return kmers.contains(kmer);
//			}else{
//				return false;
//			}
//		}catch (UnsupportedEncodingException e) {
			return this.exists(kmer);
//		}
	}
	
	public int size(){
		return kmers.size();
	}
	
	public int getKmerSize(){
		return kmerSize;
	}
	
}
