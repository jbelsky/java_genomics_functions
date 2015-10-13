package genomics_functions;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import net.sf.samtools.BAMIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;


public class BAMInput {

	public static ArrayList<String> get_chr_list(String bam_file_name){
		
		// Get the sam file
		SAMFileReader sam_file;
		
		// Check if the "*.bai" file exists
		if(new File(bam_file_name + ".bai").exists()){
				
			sam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
			
		}else{
			
			sam_file = new SAMFileReader(new File(bam_file_name));
			
		}

		List<SAMSequenceRecord> chr_list = sam_file.getFileHeader().getSequenceDictionary().getSequences();
		
		ArrayList<String> chr_str = new ArrayList<String>();
		for(int i = 0; i < chr_list.size(); i++){
			chr_str.add(chr_list.get(i).getSequenceName());
		}
		
		sam_file.close();
		
		return(chr_str);
		
	}
	
	
	public static int get_chr_length(String bam_file_name, String chr){
		
		// Get the sam file
		SAMFileReader sam_file;
		
		// Check if the "*.bai" file exists
		if(new File(bam_file_name + ".bai").exists()){
				
			sam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));
			
		}else{
			
			sam_file = new SAMFileReader(new File(bam_file_name));
			
		}
		
		// Get the header
		SAMFileHeader sam_file_header = sam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = sam_file_header.getSequence(chr).getSequenceLength();
		
		sam_file.close();
		
		return(chr_length);
		
	}
	
	public static int get_number_aligned_reads(String bam_file_name){
		
		// Get the sam file
		SAMFileReader sam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the BAMIndex
		BAMIndex bi = sam_file.getIndex();
		
		int num_align = 0;

		List<SAMSequenceRecord> chr_list = sam_file.getFileHeader().getSequenceDictionary().getSequences();
		for(int i = 0; i < chr_list.size(); i++){

			// Get the metadata for that chr
			num_align += bi.getMetaData(chr_list.get(i).getSequenceIndex()).getAlignedRecordCount();		
			
		}
		
		sam_file.close();
		
		return(num_align);
		
	}
	
	
	public static int get_number_aligned_reads(String bam_file_name, String[] chr_name){
		
		// Get the sam file
		SAMFileReader sam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the BAMIndex
		BAMIndex bi = sam_file.getIndex();
		
		int num_align = 0;
		
		for(String chr : chr_name){
		
			// Get the reference index for the chr
			int ref_idx = sam_file.getFileHeader().getSequence(chr).getSequenceIndex();
			
			// Get the metadata for that chr
			num_align += bi.getMetaData(ref_idx).getAlignedRecordCount();
		
		}
		
		sam_file.close();
		
		return(num_align);
		
	}
	
	// Returns the chromosomal length
	public static int get_number_aligned_reads_within_width(String bam_file_name, String[] chr_name,
														   	int read_width_low, int read_width_high
														   ){
		
		System.out.println("Finding the number of reads between " + read_width_low + " and " + read_width_high + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Set the read counter
		int num_reads = 0;
		
		// Iterate through each chr
		for(String chr : chr_name){

			// Open the iterator
			SAMRecordIterator bam_itr = bam_file.queryOverlapping(chr, 1, get_chr_length(bam_file_name, chr));
				
			while(bam_itr.hasNext()){
				
				// Get the reads width
				int read_width = bam_itr.next().getInferredInsertSize();
				
				// Only increment the read counter if the read's width falls within the specified range
				if(read_width >= read_width_low && read_width <= read_width_high){
					num_reads++;
				}
				
			}
		
			// Close the buffers
			bam_itr.close();
			
		}
			
		bam_file.close();
		
		System.out.println("Complete!\tThere are " + num_reads + " between " + read_width_low + " and " + read_width_high);
		
		return(num_reads);
		
	}
	
	
	
	// Returns the chromosomal length
	public static int[] get_fragment_length_distribution(String bam_file_name, 
													     int read_width_high,
													     String[] chr_dist
													    ){
		
		System.out.println("Finding the read width distribution of reads with maximum width of " + read_width_high + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Set up the storage array for the reads by width
		int[] read_width_storage = new int[read_width_high + 1];
		
		// Iterate through each chr
		for(String c : chr_dist){
			
			// Open the iterator
			SAMRecordIterator bam_itr = bam_file.queryOverlapping(c, 1, 
											get_chr_length(bam_file_name, c));
				
			while(bam_itr.hasNext()){
				
				// Get the reads width
				int read_width = bam_itr.next().getInferredInsertSize();
				
				// Only increment the read counter if the read's width falls within the specified range
				if(read_width <= read_width_high){
					read_width_storage[read_width]++;
				}
				
			}
			
			// Close the BAM itr
			bam_itr.close();
			
		}
		
		// Close the BAM file buffer
		bam_file.close();
		
		return(read_width_storage);
		
	}	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Returns the chromosomal length
	public static double[] read_in_sam_file(String bam_file_name, String chr, int shift){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Set up the storage array
		double[] chr_coverage = new double[chr_length + 1];
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		int num_reads = 0;

		int pos = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			if(read.getReadNegativeStrandFlag()){
				
				pos = read.getAlignmentEnd() - shift;
			
			}else{

				pos = read.getAlignmentStart() + shift;
				
			}
			
			// Ensure that the corrected read position falls within the chr region
			if(pos < 1){
				pos = 1;
			}else if(pos >= chr_coverage.length){
				pos = chr_coverage.length - 1;
			}
			
			// Increment the chr coverage
			chr_coverage[pos]++;
			
			// Increment the number of reads
			num_reads++;

		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	
	public static ChrStrPos[] read_in_sam_file_by_strand(String bam_file_name, String chr){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Set up the storage array
		ChrStrPos[] chr_coverage = new ChrStrPos[chr_length + 1];
		
		// Set up the new ChrStrPos
		for(int p = 1; p < chr_coverage.length; p++){
			chr_coverage[p] = new ChrStrPos();
		}
		
		// Start time
		long start_time =  System.currentTimeMillis();
				
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		int num_reads = 0;

		int pos = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			if(read.getReadNegativeStrandFlag()){
				
				// Set the pos
				pos = read.getAlignmentEnd();
			
				// Increment the chr coverage
				chr_coverage[pos].increment_neg_reads();
				
				// Increment the number of reads
				num_reads++;	
				
			}else{

				// Set the pos
				pos = read.getAlignmentStart();

				// Increment the chr coverage
				chr_coverage[pos].increment_pos_reads();
				
				// Increment the number of reads
				num_reads++;	
				
			}
						
		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	
	public static ChrStrPos[] read_in_sam_file_by_strand_density(String bam_file_name, String chr, double bw){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Create the positive and negative strand arrays
		double[] fwd_strand_arr = new double[chr_length + 1];
		double[] rev_strand_arr = new double[chr_length + 1];
		
		// Start time
		long start_time =  System.currentTimeMillis();
				
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		int num_reads = 0;

		int pos = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			if(read.getReadNegativeStrandFlag()){
				
				// Set the pos
				pos = read.getAlignmentEnd();
			
				// Set the density
				DensityEst.add_to_chr_array(rev_strand_arr, 1, pos, bw);
				
				// Increment the number of reads
				num_reads++;	
				
			}else if(!read.getReadNegativeStrandFlag()){

				// Set the pos
				pos = read.getAlignmentStart();

				// Set the density
				DensityEst.add_to_chr_array(fwd_strand_arr, 1, pos, bw);
								
				// Increment the number of reads
				num_reads++;	
				
			}
						
		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		// Set up the storage array
		ChrStrPos[] chr_coverage = new ChrStrPos[chr_length + 1];
		
		// Set up the new ChrStrPos
		for(int p = 1; p < chr_coverage.length; p++){
			chr_coverage[p] = new ChrStrPos();
			chr_coverage[p].set_pos_reads(fwd_strand_arr[p]);
			chr_coverage[p].set_neg_reads(rev_strand_arr[p]);
		}
		
		
		return(chr_coverage);
		
	}	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static int[] read_in_paired_sam_file_coverage(String bam_file_name, String chr,
														 int bp_coverage,
														 int frag_low, int frag_high
															){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Set up the storage array
		int[] chr_coverage = new int[chr_length + 1];
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		// Get the flanking region
		int flank = (bp_coverage - 1) / 2;
		
		// Set up the number of reads
		int num_reads = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			// Get the read width
			int width = read.getInferredInsertSize();
			
			// If the width is within the range, include it in the coverage map
			if((width >= frag_low) && (width <= frag_high)){
			
				// Get the start and end positions
				int read_start = read.getAlignmentStart();
				int read_end = read_start + width - 1;
			
				// Increment the coverage
				for(int p = read_start + flank; p <= read_end - flank; p++){
					chr_coverage[p]++;
				}

				// Increment the number of reads
				num_reads++;
				
			}

		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		
		// Close the BAM files
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	
	
	
	
	
	
	
	
	
	
	
	// Returns the chromosomal length
	public static HashMap<Integer, ArrayList<Integer>> 
		read_in_paired_sam_file_bed_hash(String bam_file_name, String chr,
										 int frag_low, int frag_high
										){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		// Set up the number of reads
		int num_reads = 0;
		
		// Set up the storage hash
		HashMap<Integer, ArrayList<Integer>> chr_hash = new HashMap<Integer, ArrayList<Integer>>();
		
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			// Get the read width
			int width = read.getInferredInsertSize();
			
			// If the width is within the range, include it in the coverage map
			if((width >= frag_low) && (width <= frag_high)){
			
				// Get the start and end positions
				int read_start = read.getAlignmentStart();
				
				// Check if the read_start is already in the HashMap
				if(!chr_hash.containsKey(read_start)){
					chr_hash.put(read_start, new ArrayList<Integer>());
				}
				
				// Add the read width to the list at that chromosomal position
				chr_hash.get(read_start).add(width);
				
				// Increment the number of reads
				num_reads++;
				
			}

		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		
		// Close the BAM files
		bam_itr.close();
		bam_file.close();

		return(chr_hash);
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Returns the chromosomal length
	public static BinReads[] read_in_sam_file_for_bin(String bam_file_name, String chr, 
													  int chr_start, int chr_end,
													  int bin_width, int bin_step,
													  int shift
													  ){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Set up the storage array
		BinReads[] chr_coverage = BinReads.create_bin_reads_array(chr_start, chr_end, 
																  bin_width, bin_step
																 );
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, chr_start, chr_end);
		
		int num_reads = 0;

		int pos = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			if(read.getReadNegativeStrandFlag()){
				
				pos = read.getAlignmentEnd() - shift;
			
			}else{

				pos = read.getAlignmentStart() + shift;
				
			}
			
			// Ensure that the corrected read position falls within the chr region
			if(pos < 1){
				pos = 1;
			}else if(pos > chr_length){
				pos = chr_length - 1;
			}
			
			// Increment the count over the correct BinReads
			chr_coverage[Arrays.binarySearch(chr_coverage, pos, BinReads.posCompare)].increment_count();
			
			// Increment the number of reads
			num_reads++;

		}
	
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	
	// Returns the chromosomal length
	public static BinReads[] read_in_paired_sam_file_for_bin(String bam_file_name, String chr, 
													  		 int chr_start, int chr_end,
													  		 int bin_width, int bin_step,
													  		 int frag_low, int frag_high
													  		){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Get the header
		// SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		// int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
			
		// Set up the storage array
		BinReads[] chr_coverage = BinReads.create_bin_reads_array(chr_start, chr_end, 
																  bin_width, bin_step
																 );
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, chr_start, chr_end);
		
		int num_reads = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			// Get the read coordinate
			int pos = read.getAlignmentStart();
			int width = read.getInferredInsertSize();
			
			// Adjust to the midpoint
			pos += (width/2);
			
			if(width >= frag_low && width <= frag_high){
			
				// Increment the count over the correct BinReads
				chr_coverage[Arrays.binarySearch(chr_coverage, pos, BinReads.posCompare)].increment_count();
			
				// Increment the number of reads
				num_reads++;

			}
				
		}
	
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	// Returns the chromosomal length
	public static double[] read_in_sam_file_density(String bam_file_name, String chr, 
			  										int chr_start, int chr_end,
			  										int shift, double bw
			  									    ){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Find the coverage value (normalized to 200000 reads)
		// double read_ratio = 2000000.0 / get_number_aligned_reads(bam_file_name);
		double read_ratio = 1;	
		
		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr max
		int chr_max = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Check that the chr_end isn't greater than the chr_max
		if(chr_end > chr_max)
			chr_end = chr_max;
		
		// Get the chr length
		int chr_length = chr_end - chr_start + 1;
				
		// Set up the storage array
		double[] chr_coverage = new double[chr_length];	
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryOverlapping(chr, 
															  Math.max(1, chr_start - 50),
															  Math.min(chr_end + 50, chr_max)
															 );				
		
		int num_reads = 0;

		int pos = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			if(read.getReadNegativeStrandFlag()){
				
				pos = read.getAlignmentEnd() - shift;
			
			}else{

				pos = read.getAlignmentStart() + shift;
				
			}
			
			
			// Ensure that the corrected read position falls within the chr region
			// if((pos >= (chr_start + 5000)) && (pos <= (chr_end - 5000))){
			
					if(chr.equals("12")){
						if((pos > 423475) && (pos < 500000)){
							continue;
						}
					}
				
					// Correct the read offset relative to the chr_start position
					pos -= chr_start;
					
					// Increment the chr coverage
					DensityEst.add_to_chr_array(chr_coverage, read_ratio, pos, bw);
							
					// Increment the number of reads
					num_reads++;
				
			// }
			

		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	
	
	
	
	// Returns the chromosomal length
	public static double[] read_in_paired_sam_file_density(String bam_file_name, String chr, double bw,
														   int read_width_low, int read_width_high
														  ){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Find the coverage value (normalized to 1000000 reads)
		// double read_ratio = 1000000.0 / (double) number_aligned_reads;
			
		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Set up the storage array
		double[] chr_coverage = new double[chr_length + 1];	
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		int num_reads = 0;

		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			int start = read.getAlignmentStart();
			int width = read.getInferredInsertSize();
			int pos = start + (width/2);
			
			if(width >= read_width_low && width <= read_width_high){
				
				// Increment the chr coverage
				DensityEst.add_to_chr_array(chr_coverage, 1, pos, bw);
							
				// Increment the number of reads
				num_reads++;
				
			}

		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Returns the chromosomal length
	public static double[] read_in_sam_file_chip_exo_density(String bam_file_name, String chr, double bw){
		
		System.out.println("Reading in " + bam_file_name + "...");
		
		// Get the sam file
		SAMFileReader bam_file = new SAMFileReader(new File(bam_file_name), new File(bam_file_name + ".bai"));

		// Find the coverage value (normalized to 1000000 reads)
		// double read_ratio = (1000000.0 / get_number_aligned_reads(bam_file_name)) * 2;
				
		// Get the header
		SAMFileHeader bam_file_header = bam_file.getFileHeader();
		
		// Get the chr length
		int chr_length = bam_file_header.getSequence(chr).getSequenceLength();
		
		// Set up the storage array
		double[] chr_coverage = new double[chr_length + 1];	
		
		// Start time
		long start_time =  System.currentTimeMillis();
		
		// Get the SAM Record Iterator
		// Note: Can also restrict the iterator to a particular region with ".queryContained(chr_name, start, end)"
		SAMRecordIterator bam_itr = bam_file.queryContained(chr, 1, chr_length);
		
		int num_reads = 0;

		int pos = 0;
				
		while(bam_itr.hasNext()){
		
			// Get each read (the SAMRecord object has functions to get at the SAM attributes)
			SAMRecord read = bam_itr.next();
			
			if(!read.getReadNegativeStrandFlag()){
				
				pos = read.getAlignmentStart();
				
				// Increment the chr coverage
				// DensityEst.add_to_chr_array(chr_coverage, read_ratio, pos, bw);
				chr_coverage[pos]++;	
				
				// Increment the number of reads
				num_reads++;

			}
				
		}
			
		System.out.println("\tComplete!\nRead in " + num_reads + " reads in " +
						   (System.currentTimeMillis() - start_time)/(60F * 1000F) + " minutes");
		 
		bam_itr.close();
		bam_file.close();

		return(chr_coverage);
		
	}
	
	
	
	
	
	
	
}
