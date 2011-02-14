package tools.kmer;

// hashing very much inspired from
// http://www.broadinstitute.org/annotation/conrad/api/src-html/calhoun/seq/KmerHasher.html

public class Kmer_large{

	private String kmer;
	private int mMult,lMult=5,hash=0;
	
	
	public Kmer_large(String kmer){
		this.kmer=kmer;
		mMult=(int) Math.pow(lMult, kmer.length()-1);
	}
	
	public int hashCode(){
		//optimize hashfunction for fixed length dna
		return hashFunction(kmer);
	}
	
	private int hashFunction(String s){
		if(hash==0){
			for(int i=0;i<s.length();i++){
				hash= hash*lMult + convert(s.charAt(i));
			}
		}
		return hash;
	}
	
	private int shiftHash(char nextNuc){
		return ((hash%mMult)*lMult) + convert(nextNuc);
	}
	
	public Kmer_large shift(char nextNuc){
		Kmer_large shifted= new Kmer_large(kmer.substring(1)+nextNuc);
		shifted.hash=shiftHash(nextNuc);
		return shifted;
	}
	
	private int reverseShiftHash(char nextNuc){
		return convert(nextNuc)*mMult+hash/lMult;
	}
	
	public Kmer_large reverseShift(char nextNuc){
		Kmer_large shifted= new Kmer_large(nextNuc+kmer.substring(0, kmer.length()-1));
		shifted.hash=reverseShiftHash(nextNuc);
		return shifted;
	}
	
	private int convert(char nuc){
		switch (nuc) {
		case 'A':
		case 'a':
			return 0;
		case 'T':
		case 't':
			return 1;
		case 'G':
		case 'g':
			return 2;
		case 'C':
		case 'c':
			return 3;
		case 'N':
		case 'n':
			return 4;
		default:
			System.err.println("strange nucleotide: "+nuc);
			return 6;
		}
	}
	
	public boolean equals(Object o){
		if(this==o){
			return true;
		}
		if(!(o instanceof Kmer_large)){
			return false;
		}
		Kmer_large test=(Kmer_large) o;
		return kmer.equals(test.kmer);
	}
	
	
}
