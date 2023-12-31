package omesim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

/*
 * Structure of the list of terms contributing to a trait or omics. 
 * the format of a line is composed of 10 columns (separated by "\t"):
 * 
 * #term_ID	num_cis	num_trans	trans_genes		trans_exp		traits		Weights			Model		Infinitesimal	NoiseVarianceComponent
 * ID		10		5			gene1,gene2		gene1,gene2		trt1,trt2	1.0,0.5,4.0,2.0	additive	vc_value_inf	vc_value_noise
 * 
 * Note that the the values of "trans_genes" are genes, meaning to 
 * select a few variants from that gene(s). Whereas the "trans_exp" means 
 * the expression levels of these genes are used. If the term is a gene, then
 * all these "trans_xxx" are genes from the same pathway, unless there is no 
 * pathway for the focal gene.  
 * 
 * The program does not assume that expression values of genes used to 
 * generate other genes have been assigned beforehand, and an iterative 
 * algorithm will update the missing values step-by-step.   
 * 
 * Note that the Infinitesimal will be added to the overall term without
 * considering any interactions, just like how we add the noise term.
 * 
 * The above file may be specified by the user, or generated by MainFrame.set_terms_model.
 * 
 */
public class CausalTerm {
	// basic info from the line
	public String ID; // gene ID or phenotype ID. 
	public int num_cis;
	public int num_trans;
	public String[] trans_var_genes = null; // default is null, indicating the related column in the line is "NA".
	public String[] trans_exp_genes = null;
	public String[] traits = null;
	public double[] weights;
	public String model;
	public double infinitesimal_vc_value;
	public double noise_vc_value;
	
	// data that are extracted from a MainFrame object to support the calculation  
	//
	// the first dimension is terms and the second is individuals. This is consistent to the dimensions 
	// of geno_G, exp_Z, and trait_Y in the MainFrame.java, although inconsistent to X[][] in SpecificModels.java 
	//
	// genotype locations and allele-values will be generated and extracted once and used in all iterations
	public double[][] geno_G_cis_var=null;
	public String[] geno_G_cis_loc=null;  // format: chr-index_loc_index (indexes are in MainFrame.geno_G[][]
	public double[][] geno_G_trans_var=null;
	public String[] geno_G_trans_loc=null; // format: chr-index_loc_index (indexes are in MainFrame.geno_G[][]
	// however, expressions and traits will be extracted each time as they may be changed over the iterations. 
	public double[][] exp_Z_val=null;
	public double[][] trait_Y_val=null;
	
	/*
	 * Constructor to parse the input info_line
	 * 
	 * #term_ID	num_cis	num_trans	trans_genes		trans_exp		traits		Weights			Model		Infinitesimal	NoiseVarianceComponent
	 * ID		10		5			gene1,gene2		gene1,gene2		trt1,trt2	1.0,0.5,4.0,2.0	additive	vc_value_inf	vc_value_noise
	 * 
	 * This constructor will assign genotype locations and values and, importantly, assign back to MainFrame main_frame:
	 * 		public HashMap<String, double[][]> gene2selected_var;  // map each gene to its selected expression-causal variants. 
			public HashMap<String, String[]> gene2selected_loc;    // the location of the selected expression-causal variants.
								// Format: chr-index_loc-index, the indexes are for public int[][] geno_variant_locs;
			public HashMap<String, double[][]> trait2selected_var;  // map each trait to its selected causal variants. 
			public HashMap<String, String[]> trait2selected_loc;  
	 */
	public CausalTerm(String info_line, MainFrame main_frame) {
		String[] info=info_line.split("\t");
		this.ID=info[0];
		this.num_cis=Integer.parseInt(info[1]);
		this.num_trans=Integer.parseInt(info[2]);
		if(!info[3].equals("NA")) this.trans_var_genes=info[3].split(",");
		if(!info[4].equals("NA")) this.trans_exp_genes=info[4].split(",");
		if(!info[5].equals("NA")) this.traits=info[5].split(",");
		String[] weights_string=info[6].split(",");
		this.weights=new double[weights_string.length];
		for(int t=0;t<weights_string.length;t++) 
			this.weights[t]=Double.parseDouble(weights_string[t]);
		this.model=info[7];
		 // only traits have infinitesimal_vc_value, expressions don't.
		this.infinitesimal_vc_value=isTrait(this.ID,main_frame)?Double.parseDouble(info[8]):Double.NaN;
		this.noise_vc_value=Double.parseDouble(info[9]);
		// extract genotype data from main_frame.geno_G, 
		this.sample_geno_vars_in_genes_(main_frame);
		// assign back to MainFrame
		if(isGene(this.ID, main_frame)) {
			if(this.geno_G_trans_var!=null && this.geno_G_cis_var==null) {
				main_frame.gene2selected_var.put(this.ID, this.geno_G_trans_var);
				main_frame.gene2selected_loc.put(this.ID, this.geno_G_trans_loc);
			}else if(this.geno_G_trans_var==null && this.geno_G_cis_var!=null) {
				main_frame.gene2selected_var.put(this.ID, this.geno_G_cis_var);
				main_frame.gene2selected_loc.put(this.ID, this.geno_G_cis_loc);
			}else if(this.geno_G_trans_var!=null && this.geno_G_cis_var!=null) {
				double[][] combined_var=new double[this.num_cis+this.num_trans][];
				String[] combined_loc=new String[this.num_cis+this.num_trans];
				for(int m_cis=0;m_cis<this.num_cis;m_cis++) {
					combined_var[m_cis]=this.geno_G_cis_var[m_cis];
					combined_loc[m_cis]=this.geno_G_cis_loc[m_cis];
				}
				for(int m_trans=0;m_trans<this.num_trans;m_trans++) {
					combined_var[this.num_cis+m_trans]=this.geno_G_trans_var[m_trans];
					combined_loc[this.num_cis+m_trans]=this.geno_G_trans_loc[m_trans];
				}
				main_frame.gene2selected_var.put(this.ID, combined_var);
				main_frame.gene2selected_loc.put(this.ID, combined_loc);
			}
			
		}else if(isTrait(ID, main_frame)) {
			if(this.geno_G_trans_var!=null) {
				main_frame.trait2selected_var.put(this.ID, this.geno_G_trans_var);
				main_frame.trait2selected_loc.put(this.ID, this.geno_G_trans_loc);
			}
			
		}else {
			System.out.println("Error: Term ID "+this.ID+" is not a gene or a trait!");
		}
		/* 	Note that the contributing expression and trait data will be extracted from main_frame.trait_Y and 
		 *  and main_frame.exp_Z during the iterations in calculate_this_ID(). Therefore, no operations on 
		 *  them in the constructor.
		 */  
	}

	/*
	 * load all terms from a file
	 */
	public static CausalTerm[] load_terms_from_file(String arch_detailed_file, MainFrame main_frame) {
		CausalTerm[] terms=null;
		try {
			BufferedReader br = new BufferedReader(new FileReader(arch_detailed_file));
			ArrayList<CausalTerm> terms_array=new ArrayList<CausalTerm>();
			String line=br.readLine();
			while(line.startsWith("#"))line=br.readLine();
			while(line!=null){
				terms_array.add(new CausalTerm(line, main_frame));
				line=br.readLine();
			}
			int num_terms=terms_array.size();
			terms=new CausalTerm[num_terms];
			for(int term_index=0;term_index<num_terms;term_index++) {
				terms[term_index]=terms_array.get(term_index);
			}
		}catch(Exception e) {e.printStackTrace();}
		return terms;
	}
	
	/*
	 * Based on the pre-specified causal terms, this function will:
	 * (1) collect the data and check whether they are available based on the 
	 *  main_frame.trait_finalized and main_frame.gene_exp_finalized. If all the 
	 *  terms have their flags "true", this function will assign the focal term
	 *  as "true". Note that being all zero or some values does not distinguish 
	 *  the main_frame.trait_finalized (or main_frame.gene_exp_finalized) to be "true"
	 *  or "false". This is because that the "true" value indicates that this term has
	 *  been finalized without depending on the other terms that are not finalized.
	 * (2) Form X[][] to run a model in SpecificModels using this.combine_arrays().
	 * (3) Calculate the values of this term using Model column.
	 * 
	 * Note that this function may be run many rounds to improve the number of finalized.
	 * Some terms may be never finalized, however mathematically could be converged. Genotype
	 * data will be extracted in advance (NOT in this function) therefore won't be updated 
	 * during the iterations.  
	 * 
	 * The output is "raw" values without adding infinitesimal and noise terms. 
	 * All are quantitative, not binary yet.  
	 */
	public void calculate_this_ID(MainFrame main_frame) {
		//regardless of the term being gene or trait, extract the data of contributing trans-exp and traits.
		// trans expressions and traits. gene_exp_finalized will be assigned to true if all dependent terms that are finalized
		boolean exps_finalized=true, traits_finalized=true; // default being true because that it will be true if nothing happened below (i.e., no gene-exp or traits involved) 
		if(this.trans_exp_genes!=null) {
			exps_finalized=this.get_term_vals(this.trans_exp_genes, main_frame);
		}if(this.traits!=null) {
			traits_finalized=this.get_term_vals(this.traits, main_frame);
		}
		double[][] combined_X=this.combine_arrays_weighted();
		if(isTrait(this.ID, main_frame)) {// this is a trait
			int focal_term_index=main_frame.trait_names2index.get(this.ID);
			main_frame.trait_Y[focal_term_index]=SpecificModels.response(this.model, combined_X, SpecificModels.default_weight);
			main_frame.trait_finalized[focal_term_index]= (exps_finalized && traits_finalized); // it is finalized if all depending terms are
		}else {  // must be a gene
			int focal_term_index=main_frame.gene_names2index.get(this.ID);
			main_frame.exp_Z[focal_term_index]=SpecificModels.response(this.model, combined_X, SpecificModels.default_weight);
			main_frame.gene_exp_finalized[focal_term_index]= (exps_finalized && traits_finalized); // it is finalized if all depending terms are
		}		
	}
	
	/*
	 * get genotype variants from geno_G. 
	 * 
	 * Particularly: main_frame.gene2filtered_cis_var maps gene_names to cis-variants of genes
	 * then randomly sample num_vars variants from these genes. 
	 * 
	 * Note that this method covers both cis and trans variants sampling. In case of cis,
	 * just pass only one gene (i.e., {itself}) to gene_names.  
	 * 
	 * A problem might be that the variants are NOT sampled proportional to the gene length, which is unfair for large genes. 
	 * A FUTURE development may make it more "fair" by giving more weights to larger genes.
	 */
	public void sample_geno_vars_in_genes_(MainFrame main_frame){
		if(main_frame.gene2filtered_cis_var==null) {
			System.out.println("Error: this.gene2filtered_cis_var not initiated!");
			return;
		}
		if(isGene(this.ID, main_frame)) { // this term is a gene, only which has cis
			if(this.num_cis!=-1) {
				this.geno_G_cis_var=new double[this.num_cis][];
				this.geno_G_cis_loc=new String[this.num_cis];
				double[][] vars_in_the_gene=main_frame.gene2filtered_cis_var.get(this.ID);
				String[] locs_in_the_gene=main_frame.gene2filtered_cis_loc.get(this.ID);
				for(int m_cis=0;m_cis<this.num_cis;m_cis++) {								
					int var_index=main_frame.generator.nextInt(vars_in_the_gene.length);
					this.geno_G_cis_var[m_cis]=vars_in_the_gene[var_index];
					this.geno_G_cis_loc[m_cis]=locs_in_the_gene[var_index];
				}
			}
		}
		// trans genetic variants 
		if(this.num_trans!=-1) {
			//  this.geno_G_trans_var=this.sample_geno_vars_in_genes(this.trans_var_genes, this.num_trans, main_frame);
			this.geno_G_trans_var=new double[this.num_trans][];
			this.geno_G_trans_loc=new String[this.num_trans]; 
			// distribute num_vars variants into genes in gene_names[]. the all genes get an equal share!
			for(int m_trans=0;m_trans<this.num_trans;m_trans++) {
				int gene_index=main_frame.generator.nextInt(this.trans_var_genes.length); // sample a genetic variant from this gene
				double[][] vars_in_a_gene=main_frame.gene2filtered_cis_var.get(trans_var_genes[gene_index]);
				String[] locs_in_a_gene=main_frame.gene2filtered_cis_loc.get(trans_var_genes[gene_index]);
				int var_index=main_frame.generator.nextInt(vars_in_a_gene.length);
				this.geno_G_trans_var[m_trans]=vars_in_a_gene[var_index];
				this.geno_G_trans_loc[m_trans]=locs_in_a_gene[var_index];
			}		
		}
	}
	
	/*
	 * randomly sample genotype variants from geno_G without considering genes 
	 * 
	 */
	public static double[][] sample_genotype_vars_wg(int num_sampling_vars, MainFrame main_frame){
		if(num_sampling_vars<0) // for terms without genetic component, the number of variants are set to -1 by MainFrame.set_terms_model().
			return null;
		double[][] sampled_vars=new double[num_sampling_vars][];
		for(int m=0;m<num_sampling_vars;m++) {  // no deep clone; therefore no cost of memory.
			int chr_index=main_frame.sample_chr();  // get a chr index based on the number of variants in all chrs. 
			sampled_vars[m]=main_frame.geno_G[chr_index][main_frame.generator.nextInt(main_frame.num_geno_variant_M[chr_index])];
		}
		return sampled_vars;
	}
	
	/*
	 * Getting the values of terms (gene expression or traits). 
	 * 
	 * Assign this.trait_Y_val if term_names are traits, otherwise this.exp_Z_val (term_names are genes)
	 * 
	 * Uses term_finalized[] (which will be MainFrame.gene_exp_finalized or main_frame.trait_finalized) 
	 * to check if the to-be extracted term has been finalized. If not, return a false for the focal term
	 * 
	 * Note that the function goes on even if the terms are not finalized. 
	 * It returns True if all dependent terms are finalized.  
	 * 
	 */
	public boolean get_term_vals(String[] term_names, MainFrame main_frame){
		if(!(isGene(this.ID, main_frame)||isTrait(this.ID, main_frame))) {
			System.out.println("Error: "+this.ID+" is neither a trait or a gene!");
			return false;
		}
		boolean finalized= true;
		if(isTrait(term_names[0], main_frame)) { // term_names are traits
			this.trait_Y_val=new double[term_names.length][];
			for(int t=0;t<term_names.length;t++) {
				int t_index=main_frame.trait_names2index.get(term_names[t]);
				this.trait_Y_val[t]=main_frame.trait_Y[t_index];
				if(!main_frame.trait_finalized[t_index]) finalized = false;
			}
		}else if(isGene(term_names[0], main_frame)){  //  term_names are genes
			this.exp_Z_val=new double[term_names.length][];
			for(int k=0;k<term_names.length;k++) {
				int k_index=main_frame.gene_names2index.get(term_names[k]);
				this.exp_Z_val[k]=main_frame.exp_Z[k_index];
				if(!main_frame.gene_exp_finalized[k_index]) finalized = false;
			}
		}
		return finalized;
	}
	
	/*
	 * Calculate values based on the terms for all lines.
	 * 
	 * Calculate correlation matrices BEFORE adding noise term and infinitesmal term. (the "gold-standard" relationship).
	 */
	public static void calculate_all_value(MainFrame main_frame, CausalTerm[] full_terms) {
		System.out.println("Started to iteratively calculating terms");
		// iteratively calculate the raw values of terms 
		for(int r=0;r<main_frame.iteration_rounds;r++) {
			System.out.println("===== Round "+r+" =====");
			for(int t_index=0;t_index<full_terms.length;t_index++) {
				if(isGene(full_terms[t_index].ID, main_frame) && 
						!main_frame.gene_exp_finalized[main_frame.gene_names2index.get(full_terms[t_index].ID)]) {
					full_terms[t_index].calculate_this_ID(main_frame);
					System.out.println("---- Gene "+full_terms[t_index].ID+" done ----");
				}
				if(isTrait(full_terms[t_index].ID, main_frame) && 
						!main_frame.trait_finalized[main_frame.trait_names2index.get(full_terms[t_index].ID)]) {
					full_terms[t_index].calculate_this_ID(main_frame);
					System.out.println("---- Trait "+full_terms[t_index].ID+" done ----");
				}
			}
			int genes_finalized=0, traits_finalized=0;
			for(int k=0; k<main_frame.num_gene_K; k++) genes_finalized+=main_frame.gene_exp_finalized[k]?1:0; 
			for(int t=0; t<main_frame.num_trait_T; t++) traits_finalized+=main_frame.trait_finalized[t]?1:0;
			System.out.println("Round #"+r+" Finished: "+genes_finalized+" Gene finalized; "+traits_finalized+" Traits finalized.");
		}	// Finished all rounds for biological contributions (regardless of convergence!) 
		// BEFORE adding noise or infinitesmal terms, output the correlation matrix:
		/* 	public double[][] corr_KxK;	// correlations between gene expressions 
			public double[][] corr_TxT;	// correlations between traits 
			public double[][] corr_TxK;	// correlations between traits and expressions */
		main_frame.calculate_correlations();  
		/*	public double[][] asso_MxK;	// associations between individual genetic variants and expressions 
			public double[][] asso_MxT;	// associations between individual genetic variants and traits */
		main_frame.calculate_associations(full_terms);	 
		// Next finalize them by adding noise and inf terms.
		for(int t_index=0;t_index<full_terms.length;t_index++) {
			if(isGene(full_terms[t_index].ID, main_frame)) { // adding noise term
				int focal_term_index=main_frame.gene_names2index.get(full_terms[t_index].ID);
				main_frame.exp_Z[focal_term_index]=SpecificModels.add_inf_noise_bin(
						main_frame.exp_Z[focal_term_index], main_frame, full_terms[t_index]);
			}
			if(isTrait(full_terms[t_index].ID, main_frame)){ // adding noise term and infinitesmal term (if exists); and transfer to binary (if needed).
				int focal_term_index=main_frame.trait_names2index.get(full_terms[t_index].ID);
				main_frame.trait_Y[focal_term_index]=SpecificModels.add_inf_noise_bin(
						main_frame.trait_Y[focal_term_index], main_frame, full_terms[t_index]);
			}			
		}
		System.out.println("Finished calculating terms");
	}
	
	
	
	/*
	 * Combining and rotating the data of 
	 * geno_G_cis_var, geno_G_trans_var, exp_Z_val=null, trait_Y_val,
	 * if they are not null arrays.
	 * 
	 * This will form input data, X[][], used in SpecificModels.java
	 * Note that X[][] has individuals as the first dimension, whereas the above 
	 * arrays have term first and individual next.  
	 * 
	 * TODO: double[] this.weights is multiplied in this function. 
	 */
	public double[][] combine_arrays_weighted(){		
		//sum up the number of terms and check if the sample sizes are consistent in multiple data arrays
		int num_var=0;
		int sample_size=-1;
		boolean[] initialized = new boolean[4]; // recording if the related array has been initialized
		// putting all data together for the convenience of extraction based on the initialized[]
		double[][][] all_data= new double[4][][]; 
		all_data[0]=this.geno_G_cis_var;
		all_data[1]=this.geno_G_trans_var;
		all_data[2]=this.exp_Z_val;
		all_data[3]=this.trait_Y_val;
		for(int t=0;t<4;t++) {
			if(all_data[t]!=null) {
				initialized[t]=true;
				num_var+=all_data[t].length;
				if(sample_size==-1) { // not assigned yet
					sample_size=all_data[t][0].length;
				}else { //already assigned previously
					if(sample_size!=all_data[t][0].length)
						System.out.println("Error: sample_size!=all_data[t][0].length");
				}
			}
		}
		double[][] X = new double[sample_size][num_var];
		//assign the initialized data to X[][]
		int total_term_index=0;
		for(int t=0;t<4;t++) {
			if(initialized[t]) {
				for(int k_t=0;k_t<all_data[t].length;k_t++) { // k_t is the index of terms within all_data[t]
					for(int i=0;i<sample_size;i++) {
						X[i]
								[total_term_index]=
								all_data
								[t]
										[k_t]
												[i]*
												this.weights[t];  // weights of the terms are included here
					}	
					total_term_index++;
				}
			}
		}return X;
	}
	
	public static boolean isTrait(String ID, MainFrame main_frame) {
		//return (ID.startsWith("Binary_Trait_") || ID.startsWith("Quantitative_Trait_"));
		return main_frame.trait_names2index.containsKey(ID);
	}
	
	public static boolean isGene(String ID, MainFrame main_frame) {
		return main_frame.gene_names2index.containsKey(ID);
	}
	
//	/*
//	 * write the basic info to a line
//	 */
//	public String output_line() {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
//		String output_line=this.ID+"\t";
//		output_line=output_line+this.num_cis+"\t";
//		output_line=output_line+this.num_trans+"\t";
//		for(int k=0;k<this.trans_var_genes.length;k++)
//			output_line=output_line+this.trans_var_genes[k]+((k==this.trans_var_genes.length)?"\t":",");
//		for(int k=0;k<this.trans_exp_genes.length;k++)
//			output_line=output_line+this.trans_exp_genes[k]+((k==this.trans_exp_genes.length)?"\t":",");
//		for(int t=0;t<this.traits.length;t++)
//			output_line=output_line+this.traits[t]+((t==this.traits.length)?"\t":",");
//		for(int w=0;w<this.weights.length;w++)
//			output_line=output_line+this.weights[w]+((w==this.weights.length)?"\t":",");
//		output_line=output_line+this.model+"\t";
//		output_line=output_line+this.infinitesimal_vc_value+"\t";
//		output_line=output_line+this.noise_vc_value;
//		return output_line;
//	}
	
}



