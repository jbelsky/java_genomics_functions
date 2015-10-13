package genomics_functions;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class ChIPPeakDens implements Comparable<ChIPPeakDens> {

	private String name;
	private int chr;
	private int start;
	private int end;
	private int peak_position;
	private double max_count;
	private boolean inACS = false;

	public ChIPPeakDens(String n, int c, int st, int e, int peak_pos,
			double peak_count) {
		name = n;
		chr = c;
		start = st;
		end = e;
		peak_position = peak_pos;
		max_count = peak_count;
	}

	// Return functions
	public String getName() {
		return (name);
	}

	public int getChr() {
		return (chr);
	}

	public int getStart() {
		return (start);
	}

	public int getEnd() {
		return (end);
	}

	public int getWidth() {
		return (end - start + 1);
	}

	public int getPeakPosition() {
		return (peak_position);
	}

	public double getCount() {
		return (max_count);
	}

	public boolean inACS() {
		return (inACS);
	}

	// Modifier Functions
	public void extend_peak(double c) {
		end++;
		if (c > max_count) {
			max_count = c;
			peak_position = end;
		}

	}

	public void set_acs() {
		inACS = true;
	}

	public void reverse_strand() {
		int end_temp = -end;
		end = -start;
		start = end_temp;
		peak_position = -peak_position;
	}

	public static ArrayList<ChIPPeakDens> find_subnuc_peak(double[] coverage,
			double chip_threshold, int chr, String peak_name) {

		// Create the storage ArrayList
		ArrayList<ChIPPeakDens> peak_list = new ArrayList<ChIPPeakDens>();

		// Initialize a ChIPPeak object
		ChIPPeakDens cp = new ChIPPeakDens(peak_name, 0, 0, 0, 0, 0.0);

		// Set the boolean for inPeak
		boolean inPeak = false;

		// Iterate through the coverage
		for (int p = 0; p < coverage.length; p++) {

			// Check if the value meets the threshold
			if (coverage[p] >= chip_threshold) {

				// If not in a bin, create a new chip_peak
				if (!inPeak) {
					cp = new ChIPPeakDens(peak_name, chr, p, p, p, coverage[p]);
					inPeak = true;
				} else {
					// Otherwise, extend the current ChIPPeak
					cp.extend_peak(coverage[p]);
				}

			} else if (inPeak) {
				// If the count does not meet the chip_threshold and inPeak, end
				// the peak
				peak_list.add(cp);
				inPeak = false;
			}

		}

		// If still inPeak at the end, add the final peak
		if (inPeak)
			peak_list.add(cp);

		// Ensure that there is sufficient distance between peaks
		int peak_idx = 0;

		while ((peak_idx + 1) < peak_list.size()) {

			// Get the current and subsequent peaks
			ChIPPeakDens cp_1 = peak_list.get(peak_idx);
			ChIPPeakDens cp_2 = peak_list.get(peak_idx + 1);

			if (cp_1.getEnd() + 100 > cp_2.getStart()) {

				// Find the max position and signal within both peaks
				int max_pos = cp_1.getPeakPosition();
				double max_sig = cp_1.getCount();

				if (cp_2.getCount() > max_sig) {
					max_pos = cp_2.getPeakPosition();
					max_sig = cp_2.getCount();
				}

				// Merge the peaks
				ChIPPeakDens cp_3 = new ChIPPeakDens(peak_name, cp_1.getChr(),
						cp_1.getStart(), cp_2.getEnd(), max_pos, max_sig);

				// Update the peak_list
				peak_list.set(peak_idx, cp_3);
				peak_list.remove(peak_idx + 1);

			} else {
				peak_idx++;
			}

		}

		// Ensure a sufficient peak width
		/*
		 * peak_idx = 0; while(peak_idx < peak_list.size()){
		 * if(peak_list.get(peak_idx).getWidth() < 50){
		 * peak_list.remove(peak_idx); }else{ peak_idx++; } }
		 */

		return (peak_list);

	}

	public static void write_output(BufferedWriter output,
			ArrayList<ChIPPeakDens> peak_list, int chr) throws IOException {

		// Set the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");

		// Write the output
		for (int i = 0; i < peak_list.size(); i++) {

			// Get the nucleosome
			ChIPPeakDens chip_peak = peak_list.get(i);

			output.write(chr + "," + chip_peak.getStart() + ","
					+ chip_peak.getEnd() + "," + chip_peak.getPeakPosition()
					+ "," + df.format(chip_peak.getCount()) + "\n");

		}

	}

	public static void write_output(BufferedWriter output, TF t,
			ArrayList<ChIPPeakDens> peak_list) throws IOException {

		// If the TF is on the reverse strand, invert the relative nucleosome
		// positions
		if (t.getStrand() == '-') {
			for (int i = 0; i < peak_list.size(); i++)
				peak_list.get(i).reverse_strand();
		}

		// Sort the rnp_list
		Collections.sort(peak_list);

		// Set the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");

		// Write the output
		for (int i = 0; i < peak_list.size(); i++) {

			// Get the nucleosome
			ChIPPeakDens chip_peak = peak_list.get(i);

			output.write(t.getName() + "," + t.getChr() + ","
					+ chip_peak.getStart() + "," + chip_peak.getEnd() + ","
					+ chip_peak.getPeakPosition() + ","
					+ df.format(chip_peak.getCount()) + "," + chip_peak.inACS()
					+ "\n");

		}

	}

	static Comparator<Object> posCompare = new Comparator<Object>() {

		public int compare(Object o1, Object o2) {
			ChIPPeakDens cpd = (ChIPPeakDens) o1;
			int pos = (Integer) o2;

			if (cpd.getPeakPosition() < pos) {
				return (-1);
			} else if (cpd.getPeakPosition() > pos) {
				return (1);
			} else {
				return (0);
			}

		}

	};

	@Override
	public int compareTo(ChIPPeakDens o) {
		return (this.peak_position - o.getPeakPosition());
	}

}
