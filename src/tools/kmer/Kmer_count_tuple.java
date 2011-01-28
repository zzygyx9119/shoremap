package tools.kmer;

import java.io.Serializable;

public class Kmer_count_tuple implements Comparable<Kmer_count_tuple>,Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 2L;
	private boolean[] kmer;
	private int count,hashValue;
	
	private Kmer_count_tuple(){
		hashValue=-1;
	}
	
	protected Kmer_count_tuple(boolean[] kmer, int count){
		this();
		this.kmer=new boolean[kmer.length];
		for(int i=0;i<kmer.length;i++){
			this.kmer[i]=kmer[i];
		}
		this.count=count;		
	}
	
	protected Kmer_count_tuple(String seq,int count,KmerTable_utils kt_utils){
		this();
		this.kmer=kt_utils.stringToBoolean(seq);
		this.count=count;
	}

	public boolean[] getKmer() {
		return kmer;
	}

	protected void setKmer(boolean[] kmer) {
		this.kmer = kmer;
		hashValue=-1;
	}
	
	protected int addCount(int count){
		this.count+=count;
		return getCount();
	}

	public int getCount() {
		return count;
	}

	protected void setCount(int count) {
		this.count = count;
	}
	

	public int compareTo(Kmer_count_tuple arg0) {
		if(this.kmer.length==arg0.kmer.length){
			return this.hashCode()-arg0.hashCode();
//			
//			for(int i=0;i<this.kmer.length;i++){
//				if(this.kmer[i]!=arg0.kmer[i]){
//					if(this.kmer[i]){
//						return 1;
//					}else{
//						return -1;
//					}
//				}
//			}
		}else{
			return this.kmer.length-arg0.kmer.length;
		}
//		return 0;
	}

	@Override
	public int hashCode() {
		if(hashValue==-1){
			final int prime = 2;
			int hashValue = 1;
			for(int i=0;i<kmer.length;i++){
				if(kmer[i]){
					hashValue=hashValue*prime+1;
				}else{
					hashValue*=prime;
				}
			}
			return hashValue*31+kmer.length;
		}else{
			return hashValue;
		}
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Kmer_count_tuple other = (Kmer_count_tuple) obj;
		if(this.kmer.length!=other.kmer.length)
			return false;
		if(this.hashCode()!=other.hashCode())
			return false;
		for(int i=0;i<this.kmer.length;i++){
			if(this.kmer[i]!=other.kmer[i]){
				return false;
			}
		}
		return true;
	}
	
	
}
