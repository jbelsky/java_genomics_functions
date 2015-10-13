package genomics_functions;
import java.util.HashMap;
import java.util.TreeMap;


public class ChrNameFunctions {

	public static String get_ucsc_chr_names(int chr){
		
		// Set up the chromosomal hash from integer to Roman
		TreeMap<Integer, String> chr_gb = new TreeMap<Integer, String>();
		chr_gb.put(1, "I");
		chr_gb.put(2, "II");
		chr_gb.put(3, "III");
		chr_gb.put(4, "IV");
		chr_gb.put(5, "V");
		chr_gb.put(6, "VI");
		chr_gb.put(7, "VII");
		chr_gb.put(8, "VIII");
		chr_gb.put(9, "IX");
		chr_gb.put(10, "X");
		chr_gb.put(11, "XI");
		chr_gb.put(12, "XII");
		chr_gb.put(13, "XIII");
		chr_gb.put(14, "XIV");
		chr_gb.put(15, "XV");
		chr_gb.put(16, "XVI");
		
		return(chr_gb.get(chr));
	
	}

	public static int roman_to_integer(String r){
		
		// Set up the chromosomal hash from integer to Roman
		TreeMap<String, Integer> chr_gb = new TreeMap<String, Integer>();
		chr_gb.put("I", 1);
		chr_gb.put("II", 2);
		chr_gb.put("III", 3);
		chr_gb.put("IV", 4);
		chr_gb.put("V", 5);
		chr_gb.put("VI", 6);
		chr_gb.put("VII", 7);
		chr_gb.put("VIII", 8);
		chr_gb.put("IX", 9);
		chr_gb.put("X", 10);
		chr_gb.put("XI", 11);
		chr_gb.put("XII", 12);
		chr_gb.put("XIII", 13);
		chr_gb.put("XIV", 14);
		chr_gb.put("XV", 15);
		chr_gb.put("XVI", 16);
		
		return(chr_gb.get(r));
	}
	
	public static HashMap<String, Integer> roman_to_integer_map(){
		
		// Set up the chromosomal hash from integer to Roman
		HashMap<String, Integer> chr_gb = new HashMap<String, Integer>();
		chr_gb.put("I", 1);
		chr_gb.put("II", 2);
		chr_gb.put("III", 3);
		chr_gb.put("IV", 4);
		chr_gb.put("V", 5);
		chr_gb.put("VI", 6);
		chr_gb.put("VII", 7);
		chr_gb.put("VIII", 8);
		chr_gb.put("IX", 9);
		chr_gb.put("X", 10);
		chr_gb.put("XI", 11);
		chr_gb.put("XII", 12);
		chr_gb.put("XIII", 13);
		chr_gb.put("XIV", 14);
		chr_gb.put("XV", 15);
		chr_gb.put("XVI", 16);
		
		return(chr_gb);
	}
	
	public static HashMap<String, String> sgd_numeric_to_sgd_R55(){
	
		HashMap<String, String> chr = new HashMap<String, String>();
		
		chr.put("chr01", "chr01_2006_01_20");
		chr.put("chr02", "chr02_2004_07_16");
		chr.put("chr03", "chr03_2006_01_13");
		chr.put("chr04", "chr04_2006_04_14");
		chr.put("chr05", "chr05_2000_03_16");
		chr.put("chr06", "chr06_2004_02_06");
		chr.put("chr07", "chr07_2005_12_02");
		chr.put("chr08", "chr08_2005_11_08");
		chr.put("chr09", "chr09_1994_12_10");
		chr.put("chr10", "chr10_2006_10_06");
		chr.put("chr11", "chr11_2005_12_16");
		chr.put("chr12", "chr12_2006_01_13");
		chr.put("chr13", "chr13_2004_02_27");
		chr.put("chr14", "chr14_2006_11_10");
		chr.put("chr15", "chr15_2006_01_06");
		chr.put("chr16", "chr16_2004_07_23");
		
		return(chr);
		
	}
	
}
