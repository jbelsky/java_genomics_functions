package genomics_functions;

public class ChrStrPos {
	
	double pos_reads = 0;
	double neg_reads = 0;
	
	public ChrStrPos(){
		
	}
	
	public void increment_pos_reads(){
		pos_reads++;
	}
	
	public void increment_neg_reads(){
		neg_reads++;
	}
	
	public void increment_pos_reads(double a){
		pos_reads += a;
	}
	
	public void increment_neg_reads(double a){
		neg_reads += a;
	}

	public void set_pos_reads(double a){
		pos_reads = a;
	}
	
	public void set_neg_reads(double a){
		neg_reads = a;
	}
	
	public double getPosReads(){
		return(pos_reads);
	}
	
	public double getNegReads(){
		return(neg_reads);
	}
	
	public double getTotalReads(){
		return(pos_reads + neg_reads);
	}
	
}
