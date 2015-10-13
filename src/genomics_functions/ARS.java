package genomics_functions;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

public class ARS implements Comparable<ARS> {

	protected String name;
	protected String chr;
	protected int start;
	protected int end;
			
	public ARS(String n, String c, int st, int en){
		name = n;
		chr = c;
		start = st;
		end = en;
	}
	
	// Return functions
	public String getName(){
		return(name);
	}
	
	public String getChr(){
		return(chr);
	}
		
	public int getStart(){
		return(start);
	}
	
	public int getEnd(){
		return(end);
	}
	
	// Static Functions
	public static HashMap<String, ArrayList<ARS>> read_in_tf_map(String filename, int ars_start_idx, int ars_end_idx) 
			throws IOException{
		
		// Get the Buffer for the filename and read in the header
		BufferedReader input = new BufferedReader(new FileReader(filename));
		
		// Read in the header
		String line = input.readLine();
		
		// Split the header
		String[] line_arr = line.split(",");
		
		// Set up the HashMap
		HashMap<String, ArrayList<ARS>> ars_map = new HashMap<String, ArrayList<ARS>>();
				
		while((line = input.readLine()) != null){
			
			// Split the line
			line_arr = line.split(",");
			String chr = line_arr[1];
			
			if(!ars_map.containsKey(chr)){
				ArrayList<ARS> acs_list = new ArrayList<ARS>();
				ars_map.put(chr, acs_list);
			}
			
			// Get the ARS characteristics
			String ars_name = line_arr[0];
			int start = Integer.parseInt(line_arr[ars_start_idx]);
			int end = Integer.parseInt(line_arr[ars_end_idx]);
			
			// Create a new ARS object
			ARS ars_obj = new ARS(ars_name, chr, start, end);
						
			// Enter the ARS into the map
			ars_map.get(chr).add(ars_obj);						
		
		}
		
		// Close the Buffer
		input.close();
		
		return(ars_map);
		
	}
	

	
	static Comparator<ARS> motifCompare = new Comparator<ARS>(){

		public int compare(ARS t1, ARS t2) {
			
			if(t1.getStart() < t2.getStart())
				return -1;
			else if(t1.getStart() > t2.getStart())
				return 1;
			else
				return 0;
		}
		
	};
	
	static Comparator<Object> posCompare = new Comparator<Object>(){

		public int compare(Object o1, Object o2) {
			ARS tf = (ARS) o1;
			int pos = (Integer) o2;
			
			if(tf.getStart() < pos){
				return(-1);
			}else if(tf.getStart() > pos){
				return(1);
			}else{
				return(0);
			}
			
		}
		
	};

	@Override
	public int compareTo(ARS o) {
		// TODO Auto-generated method stub
		return(this.start - o.getStart());
	}
	
	
}
