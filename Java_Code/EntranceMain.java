package omesim;

public class EntranceMain {

	public static void main(String[] args) {
		
		// user interface
		System.out.println("OmeSim: Simulating phenotype and in-between-ome(s) based on genotype");
		if(args.length==0) {
			System.out.println("Usage: ");
			System.out.println("\t-input <parameter_file>");
			System.out.println("Please check users manual for the format, structure and an example of the parameter file.");
			System.exit(0);
		}
		
		// load data
		String parameter_file=null;
		for(int p=0;p<args.length;p++) {
			if(args[p].equals("-input")) parameter_file=args[p+1];
		}
		// The MainFrame constructor below does all things:
		// assign genotype matrix this.geno_G based on filters of MAF and location_distance. 
		// read in gene information. 
		// assign cis_variants to their genes. 
		// read in pathway information 
		// setup genetic-architectures (terms and their causality)
		// load terms 
		// iteratively calculate values and generating genuine correlation arrays
		// output all files (omics values, causality graph, and genuine correlations)
		MainFrame main_frame=new MainFrame(parameter_file);
		System.out.println("ALL DONE!");

	}

}
