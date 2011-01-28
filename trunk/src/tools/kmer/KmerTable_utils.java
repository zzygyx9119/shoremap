package tools.kmer;

import java.util.BitSet;

public class KmerTable_utils {
	
	private int kmerSize;
	
	protected KmerTable_utils(int kmerSize) {
		this.kmerSize=kmerSize;
	}
	
	protected BitSet shiftBitSet(BitSet oldKmer, char nuc){
		BitSet kmer= oldKmer.get(3, kmerSize*3);
//		for(int i=0,j=3;j<kmerSize*3;i++,j++){
//			if(oldKmer.get(j)){
//				kmer.set(i);
//			}
//		}
		setPositionsInBitSet(kmer, (kmerSize-1)*3, nuc);
		kmer.set(kmerSize*3);
		return kmer;
	}
	
	protected BitSet stringToBitSet(String kmer_seq){
		BitSet kmer= new BitSet(kmerSize*3+1);
		for(int i=0,j=0;i<kmerSize;i++){
			j=setPositionsInBitSet(kmer, j, kmer_seq.charAt(i));
		}
		kmer.set(kmerSize*3);
		return kmer;
	}

	protected int setPositionsInBitSet(BitSet kmer,int pos,char nuc){
		switch (nuc) {
		case 'A':
		case 'a':
			pos+=3;
			break;
		case 'T':
		case 't':
			pos+=2;
			kmer.set(pos++);
			break;
		case 'G':
		case 'g':
			pos++;
			kmer.set(pos++);
			pos++;
			break;
		case 'C':
		case 'c':
			pos++;
			kmer.set(pos++);
			kmer.set(pos++);
			break;
		case 'N':
		case 'n':
			kmer.set(pos++);
			pos+=2;
			break;
		default:
			System.err.println("unknown char (shift) "+nuc);
			kmer.set(pos++);
			kmer.set(pos++);
			kmer.set(pos++);
			break;
		}
		return pos;
	}
	
	protected boolean[] shiftBoolean(boolean[] oldKmer,char nuc){
		boolean[] kmer= new boolean[kmerSize*3];
		//fill up remaining positions
		for(int i=0,j=3;j<kmerSize*3;i++,j++){
			kmer[i]=oldKmer[j];
		}
		setPositionsInBoolean(kmer, (kmerSize-1)*3, nuc);
		return kmer; 
	}
	
	protected boolean[] stringToBoolean(String kmer_seq){
		boolean[] kmer= new boolean[kmerSize*3];
		for(int i=0,j=0;i<kmerSize;i++){
			j=setPositionsInBoolean(kmer, j, kmer_seq.charAt(i));
		}
		return kmer;
	}
	
	protected int setPositionsInBoolean(boolean[] kmer,int pos,char nuc){
		switch (nuc) {
		case 'A':
		case 'a':
			kmer[pos++]=false;
			kmer[pos++]=false;
			kmer[pos++]=false;
			break;
		case 'T':
		case 't':
			kmer[pos++]=false;
			kmer[pos++]=false;
			kmer[pos++]=true;
			break;
		case 'G':
		case 'g':
			kmer[pos++]=false;
			kmer[pos++]=true;
			kmer[pos++]=false;
			break;
		case 'C':
		case 'c':
			kmer[pos++]=false;
			kmer[pos++]=true;
			kmer[pos++]=true;
			break;
		case 'N':
		case 'n':
			kmer[pos++]=true;
			kmer[pos++]=false;
			kmer[pos++]=false;
			break;
		default:
			System.err.println("unknown char (shift) "+nuc);
			kmer[pos++]=true;
			kmer[pos++]=true;
			kmer[pos++]=true;
			break;
		}
		return pos;
	}
}
