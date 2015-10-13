package genomics_functions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

public class Gene implements Comparable<Gene> {

	private String name;
	private String chr;
	private int start;
	private int end;
	private int tss_pos;
	private char strand;

	public Gene(String n, String c, int s, int e, char str) {
		name = n;
		chr = c;
		start = s;
		end = e;
		strand = str;
		if (strand == '+') {
			tss_pos = start;
		} else {
			tss_pos = end;
		}
	}

	// Return functions
	public String getName() {
		return (name);
	}

	public String getChr() {
		return (chr);
	}

	public int getGeneStart() {
		return (tss_pos);
	}

	public int getGeneEnd() {
		if (strand == '+') {
			return (end);
		} else {
			return (start);
		}
	}

	public char getStrand() {
		return (strand);
	}

	public int getStartCoord() {
		return (start);
	}

	public int getEndCoord() {
		return (end);
	}

	public int getGeneLength() {
		return (end - start + 1);
	}

	// Static Functions
	public static HashMap<String, ArrayList<Gene>> read_in_gene(String filename)
			throws IOException {

		// Get the Buffer for the filename and read in the header
		BufferedReader input = new BufferedReader(new FileReader(filename));

		// Set up the map
		HashMap<String, ArrayList<Gene>> gene_map = new HashMap<String, ArrayList<Gene>>();

		// Read in the header
		String line = input.readLine();

		while ((line = input.readLine()) != null) {
			String[] acs_line = line.split(",");
			String chr = acs_line[1];

			if (!gene_map.containsKey(chr)) {
				gene_map.put(chr, new ArrayList<Gene>());
			}

			// Enter in the ACS into the map
			gene_map.get(chr).add(
					new Gene(acs_line[0], chr, Integer.parseInt(acs_line[4]),
							Integer.parseInt(acs_line[5]), acs_line[3]
									.charAt(0)));

		}

		// Close the Buffer
		input.close();

		return (gene_map);

	}

	public static Comparator<Object> posCompare = new Comparator<Object>() {

		public int compare(Object o1, Object o2) {
			Gene g = (Gene) o1;
			int pos = (Integer) o2;

			if (g.getGeneStart() < pos) {
				return (-1);
			} else if (g.getGeneStart() > pos) {
				return (1);
			} else {
				return (0);
			}

		}

	};

	public static Comparator<Gene> TTSEnd = new Comparator<Gene>() {

		public int compare(Gene o1, Gene o2) {
			return (o1.getGeneEnd() - o2.getGeneEnd());
		}

	};

	@Override
	public int compareTo(Gene o) {
		// TODO Auto-generated method stub
		return (this.getGeneStart() - o.getGeneStart());
	}

}
