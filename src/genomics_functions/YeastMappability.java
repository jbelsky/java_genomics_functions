package genomics_functions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;

public class YeastMappability {

	public static void get_bin_mappability(BinReads[] bin_list, int chr) {

		// Set the mappable count header and footer
		String mappable_file_header = "c:/Users/jab112/Documents/macalpine_lab/bam_files/yeast_mapable_files/chr";
		String mappable_file_footer = "_pos_str.bam";

		// Open the relevant BAM file
		String bam_file_name = mappable_file_header + chr
				+ mappable_file_footer;
		SAMFileReader sam_file = new SAMFileReader(new File(bam_file_name));

		// Get the iterator
		SAMRecordIterator sam_itr = sam_file.iterator();

		while (sam_itr.hasNext()) {

			// Get the unmappable position
			int unmappable_pos = Integer.parseInt(sam_itr.next().getReadName());

			// Decrease the mappability in the relevant bin
			bin_list[Arrays.binarySearch(bin_list, unmappable_pos,
					BinReads.posCompare)].decrease_mappable_count();

		}

		// Close the SAM file
		sam_file.close();

		// Exclude the 10,000 bp telomeric regions from each end
		int low_telomeric_idx = Arrays.binarySearch(bin_list, 10000,
				BinReads.posCompare);
		int high_telomeric_idx = Arrays.binarySearch(bin_list,
				bin_list[bin_list.length - 1].getEnd() - 10000,
				BinReads.posCompare);

		for (int i = 0; i <= low_telomeric_idx; i++) {
			bin_list[i].set_as_unmappable();
		}

		for (int i = high_telomeric_idx; i < bin_list.length; i++) {
			bin_list[i].set_as_unmappable();
		}

		// If on chromosome 12, exclude the rDNA region
		if (chr == 12) {

			// Get the bin indices for the rDNA region
			int low_idx = Arrays.binarySearch(bin_list, 423475,
					BinReads.posCompare);
			int high_idx = Arrays.binarySearch(bin_list, 500000,
					BinReads.posCompare);

			for (int i = low_idx; i <= high_idx; i++) {
				bin_list[i].set_as_unmappable();
				System.out.println("Setting the bin beginning at "
						+ bin_list[i].getStart() + " to "
						+ bin_list[i].getEnd() + " as unmappable");
			}

		}

	}

	public static boolean[] get_mapable_regions(int chr_num, String strand_type) {

		// Get the sam file
		String mapable_bam_file_name = "c:/Users/jab112/Documents/macalpine_lab/bam_files/yeast_mapable_files/chr"
				+ chr_num + "_" + strand_type + "_str.bam";

		// Get the length for the chromosome
		int chr_length = BAMInput.get_chr_length(mapable_bam_file_name,
				Integer.toString(chr_num));

		// Create a mappable boolean
		boolean[] chr_mapable = new boolean[chr_length + 1];
		for (int k = 0; k < chr_mapable.length; k++) {
			chr_mapable[k] = true;
		}

		// Open the BAM File
		SAMFileReader mapable_sam_file = new SAMFileReader(new File(
				mapable_bam_file_name));

		// Iterate over each of the unmapable positions
		SAMRecordIterator sam_itr = mapable_sam_file.iterator();

		while (sam_itr.hasNext()) {
			chr_mapable[Integer.parseInt(sam_itr.next().getReadName())] = false;
		}

		// Close the SAM File buffer
		mapable_sam_file.close();

		// If on chromosome 12, exclude the rDNA region
		if (chr_num == 12) {
			for (int i = 423475; i <= 500000; i++) {
				chr_mapable[i] = false;
			}
		}

		// Return the boolean
		return (chr_mapable);

	}

	// Get the total mapable regions
	public static int get_total_number_mapable_regions() throws IOException {

		// Set the dir
		String dir = "c:/Users/jab112/Documents/Google Drive/Dat_Files/sacCer2_fasta/mapability/";

		// Set the total mapable number
		int total_map_num = 0;

		// Set the pattern
		Pattern p = Pattern.compile("alignment: (\\d+) \\(");

		// Iterate through each chr
		for (int i = 1; i <= 16; i++) {

			// Open the log file
			BufferedReader input = new BufferedReader(new FileReader(dir
					+ "chr" + i + "_pos_str_bowtie_" + "log_2012_09_20.txt"));

			// Read the first line
			input.readLine();

			// Get the second line
			String line = input.readLine();

			// Find the number of mapable reads
			Matcher m = p.matcher(line);

			if (m.find()) {
				total_map_num += Integer.parseInt(m.group(1));
			}

			// Close the input buffer
			input.close();

		}

		// Return the total mapable number
		return (total_map_num);

	}

}
