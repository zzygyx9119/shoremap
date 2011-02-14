package tools.kmer;

import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;

public class KmerTable {

	private ArrayList<String> files;
	private ArrayList<Kmer_count_tuple> kmers;
	private Kmer_count_tuple kmer;
	private KmerTable_utils kt_utils;
	private int curCount, limit,kmerSize;
	private String tmpPrefix;
	private ObjectOutputStream tmpOut;
	
	
	//Only to use for merging files
	public KmerTable(ArrayList<String> files){
		this(-1,"",-1);
		this.files.addAll(files);
		System.err.println("Merging...");
		
	}
	
	public KmerTable(int kmerSize,String tmpPrefix,int limit){
		this.kmerSize=kmerSize;
		this.tmpPrefix=tmpPrefix;
		this.limit=limit;
		curCount=0;
		kt_utils=new KmerTable_utils(this.kmerSize);
		files= new ArrayList<String>();
		kmers= new ArrayList<Kmer_count_tuple>();
	}
	
	public void addSequence(String seq)throws Exception{
		curCount+=1;
		kmer= new Kmer_count_tuple(seq.substring(0, kmerSize), 1, kt_utils);
		kmers.add(kmer);
		for(int i=kmerSize;i<seq.length();i++){
			kmer= new Kmer_count_tuple(kt_utils.shiftBoolean(kmer.getKmer(), seq.charAt(i)), 1);
			kmers.add(kmer);
		}
		if(curCount==limit){
			//flush
			flush();
			System.err.print("\t"+limit*files.size()+"\r");
		}
	}
	
	public void merge(String outFile) throws Exception{
		if(curCount>0){
			//write the last reads to disk
			flush();
		}
		ArrayList<ObjectInputStream> readers= new ArrayList<ObjectInputStream>();
		ArrayList<Merge_tuple> kmersOut= new ArrayList<Merge_tuple>();
		Merge_tuple mt,next_mt;
		ObjectInputStream reader;
		ObjectOutputStream out= new ObjectOutputStream(new FileOutputStream(outFile));
		for (int i=0;i<files.size();i++) {
			reader= new ObjectInputStream(new FileInputStream(files.get(i)));
			kmersOut.add(new Merge_tuple((Kmer_count_tuple) reader.readObject(), i));
			readers.add(reader);
		}
		System.err.println("merging temporary files...");
		for(int i=0;kmersOut.size()>1;i++){
			if(i%10000==0){
				System.err.print("\t"+i+"\r");
			}
			Collections.sort(kmersOut);
			mt=kmersOut.remove(0);
			kmer=mt.getKmer();
			if(mt.equals(kmersOut.get(0))){
				//Several with the same kmer ==> combine all to one
				try{
					kmersOut.add(new Merge_tuple((Kmer_count_tuple) readers.get(mt.getIndex()).readObject(),mt.getIndex()));
				}catch (EOFException e) {
					
				}
				while(kmersOut.size()>0 && mt.equals(kmersOut.get(0))){
					next_mt=kmersOut.remove(0);
					kmer.addCount(next_mt.getKmer().getCount());
					try{
						kmersOut.add(new Merge_tuple((Kmer_count_tuple) readers.get(next_mt.getIndex()).readObject(),next_mt.getIndex()));
					}catch (EOFException e) {
						
					}
				}
				out.writeObject(kmer);
			}else{
				//Only one ==> print this and continue in the file as long as the kmer is equal or smaller than
				next_mt=kmersOut.get(0);
				kmer=mt.getKmer();
				int comp=kmer.compareTo(next_mt.getKmer());
				while(comp<0){
					out.writeObject(kmer);
					try{
						kmer= (Kmer_count_tuple) readers.get(mt.getIndex()).readObject();
						comp=kmer.compareTo(next_mt.getKmer());
						if(comp>=0){
							//add it to the list
							kmersOut.add(new Merge_tuple(kmer, mt.getIndex()));
						}
					}catch (EOFException e) {
						break;
					}
				}
			}
		}
		if(kmersOut.size()>0){
			System.err.println("flushing the last reads...");
			// take the last reads from the last file
			mt=kmersOut.remove(0);
			out.writeObject(mt.getKmer());
			reader=readers.get(mt.getIndex());
			try{
				for(int i=0;true;i++){
					if(i%10000==0){
						System.err.print("\t"+i+"\r");
					}
					out.writeObject((Kmer_count_tuple)reader.readObject());
				}
			}catch (EOFException e) {
			}
		}
		out.close();
	}
	
	private void flush() throws Exception{
		//Write the current kmers to a file on the harddrive
		String curOutfile=tmpPrefix+"_kmerTable_kmers_"+(files.size()+1)+".javaObj";
		tmpOut=new ObjectOutputStream(new FileOutputStream(curOutfile));
		System.err.print("Sorting...\r");
		Collections.sort(kmers);
		System.err.print("Done sorting...\r");
		kmer=kmers.get(0);
		for(int i=0;i<kmers.size();i++){
			if(kmer.equals(kmers.get(i))){
				kmer.addCount(1);
			}else{
				//print current kmer and start with the next
				tmpOut.writeObject(kmer);
				kmer=kmers.get(i);
			}
		}
		tmpOut.writeObject(kmer);
		tmpOut.close();
		files.add(curOutfile);
		kmers=new ArrayList<Kmer_count_tuple>();
		curCount=0;
	}
}

class Merge_tuple implements Comparable<Merge_tuple>{
	
	private Kmer_count_tuple kmer;
	private int index;
	
	protected Merge_tuple(Kmer_count_tuple kmer, int index) {
		super();
		this.kmer = kmer;
		this.index = index;
	}

	protected Kmer_count_tuple getKmer() {
		return kmer;
	}

	protected int getIndex() {
		return index;
	}

	public int compareTo(Merge_tuple o) {
		return this.getKmer().compareTo(o.getKmer());
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((kmer == null) ? 0 : kmer.hashCode());
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
		Merge_tuple other = (Merge_tuple) obj;
		if (kmer == null) {
			if (other.kmer != null)
				return false;
		} else if (!kmer.equals(other.kmer))
			return false;
		return true;
	}
	
	
	
}
