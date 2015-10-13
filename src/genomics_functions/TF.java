package genomics_functions;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

public class TF implements Comparable<TF> {

	protected String name;
	protected String chr;
	protected int pos;
	protected char strand;
	protected int start;
	protected int end;
	private DecimalFormat df = new DecimalFormat("#.####");

	public TF(String n, String c, int midpoint, char str, int st, int en) {
		name = n;
		chr = c;
		pos = midpoint;
		strand = str;
		start = st;
		end = en;
	}

	// Return functions
	public String getName() {
		return (name);
	}

	public String getChr() {
		return (chr);
	}

	public int getPos() {
		return (pos);
	}

	public char getStrand() {
		return (strand);
	}

	public int getStart() {
		return (start);
	}

	public int getEnd() {
		return (end);
	}

	public int getMotifStart() {
		if (strand == '+') {
			return (start);
		} else {
			return (end);
		}
	}

	// Static Functions
	public static HashMap<String, ArrayList<TF>> read_in_tf_map(String filename)
			throws IOException {

		// Get the Buffer for the filename and read in the header
		BufferedReader input = new BufferedReader(new FileReader(filename));

		// Read in the header
		String line = input.readLine();

		// Split the header
		String[] line_arr = line.split(",");

		boolean hasStart = false;

		// Check if the 4th argument is "start"
		if (line_arr.length > 4 && line_arr[4].equals("start")) {
			hasStart = true;
		}

		// Set up the output TreeMap
		HashMap<String, ArrayList<TF>> tf_map = new HashMap<String, ArrayList<TF>>();

		while ((line = input.readLine()) != null) {
			line_arr = line.split(",");
			String chr = line_arr[1];

			if (!tf_map.containsKey(chr)) {
				ArrayList<TF> acs_list = new ArrayList<TF>();
				tf_map.put(chr, acs_list);
			}

			// Make the start and end as optional arguments
			int tf_start = 0;
			int tf_end = 0;

			if (hasStart) {
				tf_start = Integer.parseInt(line_arr[4]);
				tf_end = Integer.parseInt(line_arr[5]);
			}

			// Enter in the ACS into the map
			tf_map.get(chr).add(
					new TF(line_arr[0], chr, Integer.parseInt(line_arr[2]),
							line_arr[3].charAt(0), tf_start, tf_end));

		}

		// Close the Buffer
		input.close();

		return (tf_map);

	}

	public static ArrayList<TF> read_in_tf_list(String f) throws IOException {

		// Set up the ArrayList
		ArrayList<TF> ars_list = new ArrayList<TF>();

		// Set up the input buffer
		BufferedReader input = new BufferedReader(new FileReader(f));

		// Read in the file
		String line = input.readLine();
		String[] line_arr = line.split(",");

		boolean hasStart = false;

		// Check if the 4th argument is "start"
		if (line_arr.length > 4 && line_arr[4].equals("start")) {
			hasStart = true;
		}

		while ((line = input.readLine()) != null) {

			line_arr = line.split(",");

			// Make the start and end as optional arguments
			int tf_start = 0;
			int tf_end = 0;

			if (hasStart) {
				tf_start = Integer.parseInt(line_arr[4]);
				tf_end = Integer.parseInt(line_arr[5]);
			}

			ars_list.add(new TF(line_arr[0], line_arr[1], Integer
					.parseInt(line_arr[2]), line_arr[3].charAt(0), tf_start,
					tf_end));

		}

		input.close();
		return (ars_list);

	}

	public void write_reads_output(BufferedWriter output, double[] storage,
			int span) throws IOException {

		// Write the header info
		output.write(name + "," + chr + "," + pos + "," + strand + ",");

		String sep = ",";
		if (strand == '+') {
			for (int i = 0; i < storage.length; i += span) {
				if (i == (storage.length - 1)) {
					sep = "\n";
				}
				output.write(df.format(storage[i]) + sep);
			}
		} else {
			for (int i = storage.length - 1; i >= 0; i -= span) {
				if (i == 0) {
					sep = "\n";
				}
				output.write(df.format(storage[i]) + sep);
			}
		}

	}

	public void write_reads_output(BufferedWriter output, int[] storage,
			int span) throws IOException {

		// Write the header info
		output.write(name + "," + chr + "," + pos + "," + strand + ",");

		String sep = ",";
		if (strand == '+') {
			for (int i = 0; i < storage.length; i += span) {
				if (i == (storage.length - 1)) {
					sep = "\n";
				}
				output.write(storage[i] + sep);
			}
		} else {
			for (int i = storage.length - 1; i >= 0; i -= span) {
				if (i == 0) {
					sep = "\n";
				}
				output.write(storage[i] + sep);
			}
		}

	}

	public void write_reads_output_win(BufferedWriter output, double[] chr_vec,
			int win) throws IOException {

		// Write the header info
		output.write(name + "," + chr + "," + pos + "," + strand + ",");

		String sep = ",";
		if (strand == '+') {
			for (int i = (pos - win); i <= (pos + win); i++) {
				if (i == (pos + win)) {
					sep = "\n";
				}
				if (i >= chr_vec.length) {
					output.write("0" + sep);
				} else {
					output.write(df.format(chr_vec[i]) + sep);
				}
			}
		} else {
			for (int i = (pos + win); i >= (pos - win); i--) {
				if (i == (pos - win)) {
					sep = "\n";
				}
				output.write(df.format(chr_vec[i]) + sep);
			}
		}

	}

	public void write_reads_output_win(BufferedWriter output, int[] chr_vec,
			int win) throws IOException {

		// Write the header info
		output.write(name + "," + chr + "," + pos + "," + strand + ",");

		// Write the counts
		String sep = ",";
		if (strand == '+') {
			for (int i = (pos - win); i <= (pos + win); i++) {
				if (i == (pos + win)) {
					sep = "\n";
				}
				output.write(chr_vec[i] + sep);
			}
		} else {
			for (int i = (pos + win); i >= (pos - win); i--) {
				if (i == (pos - win)) {
					sep = "\n";
				}
				output.write(chr_vec[i] + sep);
			}
		}

	}

	public String toString() {
		String s = name + "," + chr + "," + pos + "," + strand;
		return (s);
	}

	static Comparator<TF> motifCompare = new Comparator<TF>() {

		public int compare(TF t1, TF t2) {

			if (t1.getEnd() < t2.getPos())
				return -1;
			else if (t1.getStart() > t2.getPos())
				return 1;
			else
				return 0;
		}

	};

	static Comparator<Object> posCompare = new Comparator<Object>() {

		public int compare(Object o1, Object o2) {
			TF tf = (TF) o1;
			int pos = (Integer) o2;

			if (tf.getPos() < pos) {
				return (-1);
			} else if (tf.getPos() > pos) {
				return (1);
			} else {
				return (0);
			}

		}

	};

	@Override
	public int compareTo(TF o) {
		// TODO Auto-generated method stub
		return (this.pos - o.getPos());
	}

}
