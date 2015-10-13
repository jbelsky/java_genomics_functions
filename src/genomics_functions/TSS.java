package genomics_functions;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.TreeMap;


public class TSS {

	private String name;
	private String chr;
	private int pos;
	private char strand;
		
	public TSS(String n, String c, int start, int end, char str){
		name = n;
		chr = c;
		strand = str;
		if(strand == '+'){
			pos = start;
		}else{
			pos = end;
		}
	}
	
	// Return functions
	public String getName(){
		return(name);
	}
	
	public String getChr(){
		return(chr);
	}
	
	public int getPos(){
		return(pos);
	}
	
	public char getStrand(){
		return(strand);
	}
	

	// Static Functions
	public static TreeMap<String, ArrayList<TSS>> read_in_tf(String filename) throws IOException{
		
		// Get the Buffer for the filename and read in the header
		BufferedReader acs_input = new BufferedReader(new FileReader(filename));
		acs_input.readLine();
		
		// Set up the output TreeMap
		TreeMap<String, ArrayList<TSS>> acs_map = new TreeMap<String, ArrayList<TSS>>();
		
		// Read in each ACS
		String line;
		
		while((line = acs_input.readLine()) != null){
			String[] acs_line = line.split(",");
			String chr = acs_line[1];
			
			if(!acs_map.containsKey(chr)){
				ArrayList<TSS> acs_list = new ArrayList<TSS>();
				acs_map.put(chr, acs_list);
			}
			
			// Enter in the ACS into the map
			acs_map.get(chr).add(new TSS(acs_line[0], chr,
										 Integer.parseInt(acs_line[2]),
										 Integer.parseInt(acs_line[3]),
										 acs_line[4].charAt(0)
										 ));
						
		}
		
		// Close the Buffer
		acs_input.close();
		
		return(acs_map);
		
	}
	
	public static ArrayList<TSS> read_in_tf_file(String f) throws IOException{
		
		// Set up the ArrayList
		ArrayList<TSS> ars_list = new ArrayList<TSS>();
		
		// Set up the input buffer
		BufferedReader input = new BufferedReader(new FileReader(f));
		
		// Read in the file
		String line = input.readLine();
		
		while((line = input.readLine()) != null){
		
			String[] line_arr = line.split("\t");
			ars_list.add(new TSS(line_arr[0], line_arr[1], 
								 Integer.parseInt(line_arr[2]),
								 Integer.parseInt(line_arr[3]),
								 line_arr[4].charAt(0)));
			
		}
		
		input.close();
		return(ars_list);
		
	}
	
	public void write_reads_output(BufferedWriter output, double[] chr, int win) throws IOException{
		
		DecimalFormat df = new DecimalFormat("#.####");
		
		String sep = ",";
		if(strand == '+'){
			for(int i = (pos - win); i <= (pos + win); i++){
				if(i == (pos + win)){
					sep = "\n";
				}
				output.write(df.format(chr[i]) + sep);
			}
		}else{
			for(int i = (pos + win); i >= (pos - win); i--){
				if(i == (pos - win)){
					sep = "\n";
				}
				output.write(df.format(chr[i]) + sep);
			}
		}
		
		
	}
	
	
}
