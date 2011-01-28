package tools.kmer;

import java.io.Serializable;
import java.io.UnsupportedEncodingException;
import java.util.HashSet;

import utils.BloomFilter;

public class KmerSet implements Serializable{


	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private BloomFilter<String> preFilter;
	private HashSet<String> kmers;
	private int kmerSize;
	
	public KmerSet(){
		this(-1);
	}
	
	public KmerSet(int kmerSize){
		this.kmerSize=kmerSize;
		preFilter= new BloomFilter<String>(500000000, 1000000000);
		kmers=new HashSet<String>(1000000000);
	}
	
	
	public boolean addKmer(String kmer){
		if(kmerSize==-1){
			kmerSize=kmer.length();
		}
		if(kmer.length()==kmerSize){
			try{
				preFilter.add(kmer);
			}catch (UnsupportedEncodingException e) {

			}
			return kmers.add(kmer);
		}else{
			System.err.println("Tried to add a kmer of a different length than expected:");
			System.err.println(kmer);
			return false;
		}
	}
	
	public boolean exists(String kmer){
		if(kmer.length()!=kmerSize){
			return false;
		}
		try{
			if(preFilter.contains(kmer)){
				return kmers.contains(kmer);
			}else{
				return false;
			}
		}catch (UnsupportedEncodingException e) {
			return kmers.contains(kmer);
		}
	}
	
	public int size(){
		return kmers.size();
	}
	
	public int getKmerSize(){
		return kmerSize;
	}
	
}
