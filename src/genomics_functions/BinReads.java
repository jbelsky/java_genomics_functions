package genomics_functions;

import java.util.Comparator;

public class BinReads implements Comparable<BinReads> {

	private int chr;
	private int start;
	private int end;
	private int width;
	private double count;
	private int mappable_counts;
	private boolean isARS = false;

	public BinReads(int s, int e) {
		start = s;
		end = e;
		width = end - start + 1;
		count = 0;

		// Note 2 * width for positive and negative read mappable
		mappable_counts = width;

	}

	public double getCount() {
		return (count);
	}

	public int getMid() {
		return (start + (width / 2) - 1);
	}

	public int getStart() {
		return (start);
	}

	public int getEnd() {
		return (end);
	}

	public int getWidth() {
		return (width);
	}

	public int getChr() {
		return (chr);
	}

	public double getMappability() {
		return ((double) mappable_counts / width);
	}

	public boolean isARS() {
		return (isARS);
	}

	public void reset_count() {
		count = 0;
	}

	public void increment_count() {
		count++;
	}

	public void set_count(double a) {
		count = a;
	}

	public void decrease_mappable_count() {
		mappable_counts--;
	}

	public void set_as_unmappable() {
		mappable_counts = 0;
	}

	public void set_chr(int c) {
		chr = c;
	}

	public void set_ars() {
		isARS = true;
	}

	public void update_bin_count(double[] chr_counts) {
		for (int i = start; i <= end; i++) {
			count += chr_counts[i];
		}
	}

	public double get_bin_count(double[] chr_counts) {
		double bin_counts = 0;
		for (int i = start; i <= end; i++) {
			bin_counts += chr_counts[i];
		}
		return (bin_counts);
	}

	public static BinReads[] create_bin_reads_array(int chr_length, int bin_size) {

		// Get the optimal number of bins
		int num_bins;

		if (chr_length % bin_size != 0) {
			num_bins = (chr_length / bin_size) + 1;
		} else {
			num_bins = (chr_length / bin_size);
		}

		// Set up the BinReads array
		BinReads[] storage = new BinReads[num_bins];

		// Iterate through each element of the array and set up the storage
		int pos = 1;
		for (int i = 0; i < storage.length; i++) {
			storage[i] = new BinReads(pos, pos += (bin_size - 1));
			pos++;
		}

		// Set the end of the final BinReads as the chr_length
		storage[storage.length - 1].end = chr_length;

		return (storage);

	}

	public static BinReads[] create_bin_reads_array(int chr_length,
			int bin_size, int step_size) {

		// Get the optimal number of bins
		int num_bins;

		if (chr_length % step_size != 0) {
			num_bins = (chr_length / step_size) + 1;
		} else {
			num_bins = (chr_length / step_size);
		}

		// Set up the BinReads array
		BinReads[] storage = new BinReads[num_bins];

		// Iterate through each element of the array and set up the storage
		int pos = 1;
		for (int i = 0; i < storage.length; i++) {
			storage[i] = new BinReads(pos, Math.min(pos + (bin_size - 1),
					chr_length));
			pos += step_size;
		}

		// Set the end of the final BinReads as the chr_length
		storage[storage.length - 1].end = chr_length;

		return (storage);

	}

	public static BinReads[] create_bin_reads_array(int chr_start, int chr_end,
			int bin_size, int step_size) {

		// Get the width of the region of interest
		int roi_width = chr_end - chr_start + 1;

		// Get the optimal number of bins
		int num_bins;

		if (roi_width % step_size != 0) {
			num_bins = (roi_width / step_size) + 1;
		} else {
			num_bins = (roi_width / step_size);
		}

		// Set up the BinReads array
		BinReads[] storage = new BinReads[num_bins];

		// Iterate through each element of the array and set up the storage
		int pos = chr_start;
		for (int i = 0; i < storage.length; i++) {
			storage[i] = new BinReads(pos, Math.min(pos + (bin_size - 1),
					chr_end));
			pos += step_size;
		}

		return (storage);

	}

	public static BinReads[] create_bin_reads_array_around_feature_position(
			int pos, int feature_width, int bin_width) {

		// Get the number of bins
		int num_bins = 2 * feature_width / bin_width;

		// Set up the storage array
		BinReads[] storage = new BinReads[num_bins];

		// Set the starting position
		int bin_start = pos - feature_width;

		// Create the positions for the bins
		for (int i = 0; i < storage.length; i++) {

			// Get the bin_end position
			int bin_end = bin_start + bin_width - 1;

			// Create the bin
			storage[i] = new BinReads(bin_start, bin_end);

			// Update the bin_start
			bin_start += bin_width;

			// Remove the bin around the feature position
			if (bin_start == pos) {
				bin_start++;
			}

		}

		// Return the array
		return (storage);

	}

	public static BinReads[] create_bin_reads_dist(int chr_start, int chr_end,
			int num_bins) {

		// Get the width of the region of interest
		int roi_width = chr_end - chr_start + 1;

		// Get the step size
		int step_size = roi_width / num_bins;

		// Get the left_over
		int remaining_dist = roi_width % num_bins;

		// Set up the BinReads array
		BinReads[] storage = new BinReads[num_bins];

		// Iterate through each element of the array and set up the storage
		int pos = chr_start;
		for (int i = 0; i < storage.length; i++) {

			// Get the step size
			int bin_step = step_size;

			// Get the step size
			if (remaining_dist > 0) {
				bin_step++;
				remaining_dist--;
			}

			storage[i] = new BinReads(pos, pos + bin_step - 1);
			pos += bin_step;
		}

		return (storage);

	}

	public static Comparator<Object> posCompare = new Comparator<Object>() {

		public int compare(Object o1, Object o2) {
			BinReads br = (BinReads) o1;
			int pos = (Integer) o2;

			if (br.getEnd() < pos) {
				return (-1);
			} else if (br.getStart() > pos) {
				return (1);
			} else {
				return (0);
			}

		}

	};

	public int compareTo(BinReads br_compare) {

		// If the end of this BinReads is less than the start of br_compare,
		// this bin comes first
		if (this.end < br_compare.getStart()) {
			return (-1);
		} else if (this.start > br_compare.getEnd()) {
			// If the start of this BinReads is greater than the end of
			// br_compare,
			// this bin comes later
			return (1);
		} else {
			// Otherwise, there is overlap between the bins
			return (0);
		}

	}

}
