package tools.fasta;

import java.io.FileNotFoundException;
import java.io.IOException;

import rfBio.*;

public class split_multifasta {

	public static void main(String[] args)
	 throws IOException, FileNotFoundException{
		

//		Bioutil test = new Bioutil();
		
		if(args.length != 3){
			System.out.println("Wrong number of parameters!" + '\n' + "Usage: split_multifasta infile outfile num_seq_per_file.");
			System.exit(0);
			
		}
		fasta f = new fasta();
		f.splitLarge(args[0], args[1], Type.toInt(args[2]));
	}

}
