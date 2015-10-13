package genomics_functions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class PWMFunctions {

	@SuppressWarnings("serial")
	private static final HashMap<Character, Integer> nucleotide_to_index = new HashMap<Character, Integer>() {
		{
			put('A', 0);
			put('C', 1);
			put('G', 2);
			put('T', 3);
		}
	};

	public static HashMap<Character, Integer> get_nuc_to_arr_index_hash() {
		return (nucleotide_to_index);
	}

	public static HashMap<String, Double> read_in_bg_model(
			String bg_model_file_name) throws IOException {

		// Create the HashMap
		HashMap<String, Double> bg_map = new HashMap<String, Double>();

		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(
				bg_model_file_name));

		String line;

		while ((line = input.readLine()) != null) {
			String[] line_arr = line.split("\t");
			String kmer = line_arr[0];
			if (!kmer.equals("kmer")) {
				double score = Double.parseDouble(line_arr[1]);
				bg_map.put(kmer, score);
			}
		}

		input.close();

		return (bg_map);

	}

	// Read in the PWM file
	public static double[][] read_in_pwm(String pwm_file) throws IOException {

		// Create the double arraylist
		ArrayList<Double[]> pwm = new ArrayList<Double[]>();

		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(pwm_file));

		// Read in the header
		String line = input.readLine();

		while ((line = input.readLine()) != null) {

			// Split the line on tab
			String[] line_arr = line.split("\t");

			// Create the double array
			Double[] pwm_arr = new Double[4];

			for (int i = 0; i < 4; i++) {
				pwm_arr[i] = Double.parseDouble(line_arr[i]);
			}

			// Add this double array to the ArrayList
			pwm.add(pwm_arr);

		}

		// Create the double array from the arraylist
		double[][] pwm_dbl_arr = new double[4][pwm.size()];
		for (int i = 0; i < pwm.size(); i++) {
			Double[] cur_pos_arr = pwm.get(i);
			for (int j = 0; j < 4; j++) {
				pwm_dbl_arr[j][i] = cur_pos_arr[j];
			}
		}

		input.close();

		return (pwm_dbl_arr);

	}

	// Read in the PWM file
	public static Double[][] read_in_transfac_pwm(String pwm_file)
			throws IOException {

		// Create the double arraylist
		ArrayList<Double[]> pwm = new ArrayList<Double[]>();

		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(pwm_file));

		String line;

		while ((line = input.readLine()) != null) {

			// Split the line on tab
			String[] line_arr = line.split("\t");

			// Create the double array
			Double[] pwm_arr = new Double[line_arr.length];

			for (int i = 0; i < pwm_arr.length; i++) {
				pwm_arr[i] = Double.parseDouble(line_arr[i]);
			}

			// Add this double array to the ArrayList
			pwm.add(pwm_arr);

		}

		// Create the double array from the arraylist
		Double[][] pwm_dbl_arr = new Double[4][pwm.get(0).length];
		for (int i = 0; i < pwm.size(); i++) {
			Double[] cur_pos_arr = pwm.get(i);
			for (int j = 0; j < cur_pos_arr.length; j++) {
				pwm_dbl_arr[i][j] = cur_pos_arr[j];
			}
		}

		input.close();

		return (pwm_dbl_arr);

	}

	public static double[][] reverse_complement_pwm(double[][] pwm) {

		// Set up the storage pwm
		double[][] pwm_rc = new double[4][pwm[0].length];

		// Update the pwm_rc
		for (int i = 0; i < pwm[0].length; i++) {

			for (int j = 0; j < 4; j++) {

				// Update the reverse complement pwm
				pwm_rc[j][i] = pwm[(3 - j)][(pwm[0].length - 1 - i)];

			}

		}

		// Return the pwm_rc
		return (pwm_rc);

	}

	public static double find_bg_score(String seq, int width,
			HashMap<String, Double> bg_score_map) {

		// Initialize the parameters
		double score = 0;
		int lower = 0;
		int w = 1;

		// Iterate through the sequence to find the score
		for (int i = 0; i < seq.length(); i++) {

			// Increment the log-likelihood score
			score += Math.log10(bg_score_map.get(seq
					.substring(lower, lower + w)));

			// If the current width is below the bg order width, increment the
			// width
			if (w < 4) {
				w++;
			} else {
				// Once the current width has reached the bg order width,
				// increment the position
				lower++;
			}
		}

		return (score / Math.log10(2));

	}

	public static double find_pwm_score(String seq, double[][] pwm) {

		// Intialize the score
		double score = 0;

		for (int p = 0; p < seq.length(); p++) {

			// Get the probability at the position
			double position_prob = pwm[nucleotide_to_index.get(seq.charAt(p))][p];

			// Set the position probability if 0 to 0.001
			if (position_prob < 0.001)
				position_prob = 0.001;

			// Update the score
			score += Math.log10(position_prob);

		}

		return (score / Math.log10(2));
	}

	public static double get_log_score(String chr_seq, double[][] pwm, int p,
			char str, HashMap<String, Double> bg_model) {

		// Get the seq
		String seq = get_seq_to_evaluate(chr_seq, p, pwm[0].length, str);

		// Calculate the background score
		double bg_score = find_bg_score(seq, 4, bg_model);

		// Calculate the ACS score
		double acs_score = find_pwm_score(seq, pwm);

		// Calculate the log-odds score
		double log_diff_score = acs_score - bg_score;

		// Return the log_diff_score
		return (log_diff_score);

	}

	public static String get_seq_to_evaluate(String chr_seq, int pos,
			int pwm_width, char str) {

		// Set the seq variable
		String seq = "";

		// Get the seq on either the forward or reverse strand
		if (str == '+') {
			seq = chr_seq.substring(pos, pos + pwm_width);
		} else {
			seq = GenomeFunctions.reverse_complement_sequence(chr_seq
					.substring(pos - (pwm_width - 1), pos + 1));
		}

		// Return the sequence
		return (seq);

	}

}
