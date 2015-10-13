package genomics_functions;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class GenomeFunctions {
	
	@SuppressWarnings("serial")
	private static final HashMap<Character, Character> nuc_comp = 
		new HashMap<Character, Character>(){{
		put('A', 'T');
		put('C', 'G');
		put('G', 'C');
		put('T', 'A');
		put('N', 'N');
	}};
		
	public static String read_in_chr_fasta_file(String c, String chr_header) throws IOException{
		
		// Read in the chromosome fasta file
		BufferedReader fasta_file = new BufferedReader(new FileReader(chr_header + c + ".fasta"));
		
		// Read the header
		fasta_file.readLine();
		
		// Read in the chr and make it a StringBuffer
		StringBuffer chr = new StringBuffer();
		
		// Add in an "N" at position 0 to offset the buffer
		chr.append("N");
		
		String line;
		System.out.println("Reading in the chr " + c);
		while((line = fasta_file.readLine()) != null){
			chr.append(line);
		}
		fasta_file.close();
		
		System.out.println("\tComplete!");		
		
		return(chr.toString());
		
	}
	
	public static String read_in_chr_split_fasta_file(String c, String chr_header, int file_idx) throws IOException{
		
		// Create the file footer
		String file_footer = "chr" + c + "/chr" + c + "_" + file_idx + ".fasta";
		
		// Read in the chromosome fasta file
		BufferedReader fasta_file = new BufferedReader(new FileReader(chr_header + file_footer));
		
		// Read the header
		fasta_file.readLine();
		
		// Read in the chr and make it a StringBuffer
		StringBuffer chr = new StringBuffer();
		
		// Add in an "N" at position 0 to offset the buffer
		chr.append("N");
		
		String line;
		System.out.println("Reading in the chr " + c);
		while((line = fasta_file.readLine()) != null){
			chr.append(line);
		}
		fasta_file.close();
		
		System.out.println("\tComplete!");		
		
		return(chr.toString());
		
	}
	

	public static String reverse_complement_sequence(String seq){
        
        // Reverse the String
        String seq_rev = new StringBuffer(seq).reverse().toString();

        // Create the output StringBuilder
        StringBuilder seq_rc = new StringBuilder();
        for(int j = 0; j < seq_rev.length(); j++){
                seq_rc.append(nuc_comp.get(seq_rev.charAt(j)));
        }   

        return(seq_rc.toString());

    }
	
	
}
