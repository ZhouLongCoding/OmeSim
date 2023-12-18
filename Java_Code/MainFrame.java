package omesim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

/*
 * Main Constructor parsing all the parameters:
 * 1. Genotype data and QC parameters
 * 2. High-level architectures
 * 3. Specific statistic models
 * 4. Output folder and files
 *
 * The parameters are specified in the format of XYZ=values. 
 */
public class MainFrame {
	
	public static final int random_seed =1;
	public final Random generator = new Random(MainFrame.random_seed);
	
	public static final String[] supported_genotype_file_format= {"CSV","VCF","PLINK"};
	// files' path info
	public String genotype_file="";
	public String genotype_file_format="";
	public String genes_file=""; 
	public String pathway_file="";
	public String output_folder="";
	
	// pop-gen parameters 
	public double max_maf=0.5;
	public double min_maf=0.05;
	public int location_distance=0; // this is used to guard against LD;
	
	// genetic architecture and model input parameters (all parameters are general proportions to guide the random assignment of terms)
	// based on the probabilities of contributors below, each gene/trait will be randomly decided on whether it contains corresponding contributors 
	public double[] gene_contributors= {1.0, 0.5, 0.3, 0};  // proportion of gene expressions that have the contribution from [0]=cis;[1]=trans;[2]=other_gene_exp;[3]=trait.  
	public double[] traits_contributors={1.0, 0.8, 0.5, 0.3};  //proportion of traits that have the contribution from [0]=genetics;[1]=gene_exp;[2]=trait,[3]=infinitesimal.
	// parameters to guide percentage of contributions 
	public double[] gene_contri_weights= {1.0, 0.5, 4.0, 2.0};  // weights of gene-contributors from [0]=cis;[1]=trans;[2]=other_gene_exp;[3]=trait.  
	public double[] traits_contri_weights={1.0, 4.0, 2.0};  // weights of trait-contributors from [0]=genetics;[1]=gene_exp;[2]=trait; No infinitesimal at this stage as [4]=infinitesimal will be added later before adding noise.
	public double weights_relative_range = 0.5; 			// relative range of the above gene_contri_weights[] and traits_contri_weights[]
	public double traits_var_comp_inf_min = 0.05;  // minimal variance component for the infinitesimal term: traits only (no genes).  
	public double traits_var_comp_inf_max = 0.45;  // maximal variance component for the infinitesimal term: traits only (no genes).  
	public double var_comp_noise_min = 0.1;  // Note that the non-noise part may not be deemed as "heritability" as the other traits (and noises in the traits and expressions) are also involved.  
	public double var_comp_noise_max = 0.9;  // Note that the noise term will be added after the infinitesimal term, and the nomalization will be conducted again.   

	public int cis_var_flanking = 500000;	// flanking regions in base-pair defining the region for cis-variants. 
	public int cis_variant_numbers_mean= 5;  // mean of number of cis-variants of a gene 
	public int cis_variant_numbers_range= 10;  // range of number of cis-variants of a gene. Note that negative number will be converted to zero 
	public int trans_variant_numbers_mean= 5;  // mean of number of trans-variants of a gene (multiple genes will share the total number of variants.)
	public int trans_variant_numbers_range= 10;  // range of number of trans-variants of a gene (multiple genes will share the total number of variants.)
	public int num_contributing_genes_mean=4; 	 // mean of number of trans-genes contributing to exp or traits. Note that the SAME parameter controls both genetic and exp. 
	public int num_contributing_genes_range=4; 	 // range of number of trans-genes contributing to exp or traits. Note that the SAME parameter controls both genetic and exp. 
	public int num_contributing_traits_mean=2; 	 // mean of number of traits contributing to exp or traits. 
	public int num_contributing_traits_range=1; 	 // range of number of traits contributing to exp or traits. 
	public int num_infinitesimal = 5000; 		// the number of genetic variants whose linearly added effect will serve as the infinitesimal term 
	// based on the probabilities below, each gene/trait will be randomly decided on which model it follows. (The sum of all probabilities must be 1.0).  
	public double[] gene_models= {0.5,0.05,0.05,0.1,0.3}; // proportion of gene expressions with models [0]=additive;[1]=epistatic;[2]=compensatory; [3]=heterogenous; [4]=compound 		
	public double[] trait_models= {0.5,0.05,0.05,0.1,0.3}; // proportion of traits with models [0]=additive;[1]=epistatic;[2]=compensatory; [3]=heterogenous; [4]=compound
	
	public int iteration_rounds= 10;  // number of iterations assigning all terms. (Iteration is needed because of that many terms are not initialized. 
	
	public double trait_binary_proportion=0.5;   // proportion of traits that are binary. The rest are quantitative.
	public static double[] binary_mode_proportion= {0.5,0.5};  // proportion of modes to decide binary: [0]=liability; [1]=logistic
	public static String[] binary_modes= {"liability", "logistic"};

	// Data structures that will be generated based on the input files 
	public int num_trait_T=10;		// #traits T default = 10
	public String[] trait_names;// in the format of T_ID;
	public HashMap<String, Integer> trait_names2index = new HashMap<String, Integer>();
	public boolean[] trait_finalized;  // finalized means no dependence on other terms that are not finalized. 
		//Initially, in the first round of calculation, terms only relying on genetics will be finalized. Iteratively, many others will be finalized too. 
	
	public int num_gene_K;		// #genes K
	public String[] gene_names;	// just copied from the genes_file
	public HashMap<String, Integer> gene_names2index = new HashMap<String, Integer>();
	public HashMap<String, int[]> gene_locs = new HashMap<String, int[]>(); //chr, start and end locations of genes
	public boolean[] gene_exp_finalized; // finalized means no dependence on other terms that are not finalized. 
		//Initially, in the first round of calculation, terms only relying on genetics will be finalized. Iteratively, many others will be finalized too. 
	
	public int num_subj_N;		// #subjects, or sample size, N. 
	public String[] subj_names;	// names of all subjects
	int num_chr;				// number of chromosomes
	String[] chr_names;			// names of chromosomes
	HashMap<String, Integer> chr2index=new HashMap<String, Integer>(); // chr names to index (starts from Zero) 
	public int[] num_geno_variant_M;	// #SNPs M for all chromosomes
	public int[][] geno_variant_locs;	// the locations of the above num_geno_variant_M
	public double[][][] geno_G;	// genotype array: chr x loc x individual 
								// may be generated by randomly selecting variants from the genotype file
	public HashMap<String, double[][]> gene2filtered_cis_var;  // map each gene to its cis-variants within the flanking region filtered by MAF and distances (LD).
	
	public String[] pathway_names;
	public HashMap<String, String[]> pathways2genes = new HashMap<String, String[]>(); // pathway => its members
	public HashMap<String, ArrayList<String>> gene2pathways = new HashMap<String, ArrayList<String>>(); // gene => its pathways
	 
	//public double[][] GRM;	// GRM between individuals based on genetics
								// for the function SpecificModels.add_infinitesimal_term()
								// This may be generated by the full genome, instead of the selected geno_G above. 
								// NOTE: It is not supported in the first version (2023-July). Currently using many trans-variants to represent it.
	
	// genetic architecture file generated based on the genetic architecture input parameters
	String arch_detailed_file;
	
	// Data structures for output. These are gold-standard for verifying statistic models/tools discovering them
	// Actual exp and trait generated
	public double[][] exp_Z; 	// exp x individual; could be gene expressions or other omics data (e.g., protein)
	public double[][] trait_Y; 	// trait x individual; could be phenotypic traits 
	// variance components contributing to exp and trait from different categories. 
	// Note that these are based on a linear regression between the focal trait/exp and the category terms. It is an approximation. Also, the sum of all categories may not be 1.0.
	// The above regression is still considered part of "gold-standard" because the focal trait/exp are before adding the noise term.  
	public double[][] var_comp_exp=new double[4][];   // variance components of expressions:[0]=cis;[1]=trans;[2]=other_gene_exp;[3]=trait
	public double[][] var_comp_trait=new double[3][]; //variance components of traits:[0]=genetics;[1]=gene_exp;[2]=trait.
	// correlation structures before adding the noise terms. 
	public double[][] corr_KxK;	// correlations between gene expressions, a triangle matrix this.corr_KxK[k1][k2-k1] = corr(k1,k2) 
	public double[][] corr_TxT;	// correlations between traits, a triangle matrix  this.corr_TxT[t1][t2-t1] = corr(t1,t2) 
	public double[][] corr_TxK;	// correlations between traits and expressions
	public double[][] asso_MxK;	// associations between individual genetic variants and expressions 
	public double[][] asso_MxT;	// associations between individual genetic variants and traits
								// The "association" is also actually correlation. We name it differently to emphasize it is between genetics and terms  
	
	public HashMap<String, ArrayList<String>> causality_graph;  // a directed graph recording who generates whom;
								// each key is a node, and the value is an ArrayList of the terms contribute to the node.  		
	
	public MainFrame(String parameter_file) {
		try {
			BufferedReader br= new BufferedReader(new FileReader(parameter_file));
			String line=br.readLine();
			while(line!=null) {
				// skip headers
				if(line.startsWith("#"))line=br.readLine(); 
				// parse parameters from the file
				String[] para=line.split("="); //para[0]=name; para[1]=value.
				// extract file paths of real-data  				
				if(para[0].equals("Genotype_File"))this.genotype_file=para[1];
				if(para[0].equals("Genotype_File_Format"))this.genotype_file_format=para[1];
				if(para[0].equals("Genes_File"))this.genes_file=para[1];
				if(para[0].equals("Pathway_File"))this.pathway_file=para[1];
				// extract file paths of detailed genetic architecture generated by this class. 
				if(para[0].equals("Arch_Detailed_File"))this.arch_detailed_file=para[1];
				// how many traits to be simulated [default = 10]
				if(para[0].equals("Num_Traits")) this.num_trait_T=Integer.parseInt(para[1]);
				// extract pop-gen parameters
				if(para[0].equals("Max_MAF"))this.max_maf=Double.parseDouble(para[1]);
				if(para[0].equals("Min_MAF"))this.min_maf=Double.parseDouble(para[1]);
				if(para[0].equals("Loc_Distance"))this.location_distance=Integer.parseInt(para[1]);
				// extract detailed models
				if(para[0].equals("Genotype"))this.genotype_file=para[1];
				if(para[0].equals("Genotype"))this.genotype_file=para[1];
				if(para[0].equals("Genotype"))this.genotype_file=para[1];
				// output folder
				if(para[0].equals("Output_Folder"))this.output_folder=para[1]+"/";

				line=br.readLine();
			}br.close();
			// assign genotype matrix this.geno_G based on filters of MAF and location_distance. 
			// supporting three file types VCF, CSV, and PLINK (tped)
			this.readin_genotype_with_maf_ld_filters();
			// read in gene information. 
			// Note that it relies on the this.readin_genotype_with_maf_ld_filters() to set up chr2index 
			this.read_in_genes_locs();
			// assign cis_variants to their genes. 
			this.set_gene2cis_var_map();
			// read in pathway information 
			this.read_in_pathways();
			// setup genetic-architectures 
			this.set_terms_model();
			// load terms generated above
			CausalTerms[] arch_terms=CausalTerms.load_terms_from_file(this.arch_detailed_file); 
			// iteratively calculate
			CausalTerms.calculate_all_value(this, arch_terms);			
			// output all files (omics values, causality graph, and genuine correlations) 
			this.output_files(this.output_folder);
			
		}catch(Exception e) {e.printStackTrace();}
	}
	
	/*
	 * assign genotype matrix this.geno_G based on filters of MAF and location_distance. 
	 * supporting three file types VCF, CSV, and PLINK (tped)
	 * 
	 * It will initiate:
	 * 		public int num_subj_N;		// #subjects, or sample size, N. 
			public String[] subj_names;	// names of all subjects
			int num_chr;				// number of chromosomes
			HashMap<String, Integer> chr2index=new HashMap<String, Integer>(); // chr names to index (starts from Zero) 
			String[] chr_names;			// name of chromosomes
			
			public int[] num_geno_variant_M;	// #SNPs M for all chromosomes
			public int[][] geno_variant_locs;	// the locations of the above num_geno_variant_M
			public double[][][] geno_G;	// genotype array: chr x loc x individual
	 * 
	 * It filter variants based on their MAF the distances between them. 
	 * 
	 * The CSV file should be formatted as: 
	 * ## infromation headers
	 * #chr,loc,subj1,subj2,subj3,...,subjN
	 * 1,1000,1,2,0,...,1.
	 * 
	 * Note that the alleles is coded as 012, instead of ATCG. 
	 * Other numerical codings (e.g., incorporating rep-learn scores) are also supported;
	 * however, the definition of MAF might be questionable.  
	 */
	public void readin_genotype_with_maf_ld_filters() {
		try {
			BufferedReader br=new BufferedReader(new FileReader(this.genotype_file));
			String line=br.readLine();
			ArrayList<String> chr_names_array=new ArrayList<String>();
			ArrayList<Integer> num_geno_variant_array=new ArrayList<Integer>();			
			ArrayList<ArrayList<Integer>> geno_variant_locs_array=new ArrayList<ArrayList<Integer>>(); 
			ArrayList<ArrayList<double[]>> geno_G_array=new ArrayList<ArrayList<double[]>>(); 
			int last_loc=-1;		// this -1 is only an indicator that the last_loc hasn't been assigned. 
			int current_chr_index=-1;  // started with -1, and will be updated at the first line to 0. 
			String current_chr_name="PlaceHolderThatIsNotChrName";  
			if(this.genotype_file_format.equals("CSV")) {
				while(line.startsWith("##")) { // skip the general headers, if any
					line=br.readLine();
				}
				String[] header=line.split(",");
				this.num_subj_N=header.length-2;
				this.subj_names=new String[this.num_subj_N];
				for(int n=0;n<this.num_subj_N;n++)this.subj_names[n]=header[n+2];
				while(line!=null) {
					String[] tmp=line.split(",");
					if(!current_chr_name.equals(tmp[0])) { // a new chromosome					
						last_loc=0; 
						current_chr_index++;
						current_chr_name=tmp[0];
						this.chr2index.put(current_chr_name, current_chr_index);
						chr_names_array.add(current_chr_name);
						//move forward until identifying the first line with a qualified MAF
						double[] alleles=new double[this.num_subj_N];
						for(int n=0;n<this.num_subj_N;n++)alleles[n]=Double.parseDouble(tmp[n+2]);
						double maf= maf(alleles);
						while (!(maf>=this.min_maf && maf <=this.max_maf)){
							line=br.readLine();
							tmp=line.split(",");
							if(!tmp[0].equals(current_chr_name)) {
								System.out.println("Error: No genetic variant has qualified MAF in Chr "+current_chr_name+"!");
								return;
							}
							for(int n=0;n<this.num_subj_N;n++)alleles[n]=Double.parseDouble(tmp[n+2]);
							maf= maf(alleles);	
						} // Assuming the above iteration will end (i.e., every chromosome has at least one variant with qualified MAF) 
						last_loc = Integer.parseInt(tmp[1]);
						num_geno_variant_array.add(1);  // number of variants per chr
						ArrayList<Integer> this_chr_locs=new ArrayList<Integer>();  // locations of variants per chr
						this_chr_locs.add(Integer.parseInt(tmp[1]));
						geno_variant_locs_array.add(this_chr_locs);
						ArrayList<double[]> this_chr_geno_G=new ArrayList<double[]>(); // genotype arrays per chr
						this_chr_geno_G.add(alleles);
						geno_G_array.add(this_chr_geno_G);									
					}else {		// the same chromosome, continue adding variants
						int current_loc = Integer.parseInt(tmp[1]);
						if(current_loc - last_loc>=this.location_distance) { // satisfied the distance (LD) interval.
							double[] alleles=new double[this.num_subj_N];
							for(int n=0;n<this.num_subj_N;n++)alleles[n]=Double.parseDouble(tmp[n+2]);
							double maf= maf(alleles);
							if(maf>=this.min_maf && maf <=this.max_maf) { // satisfied the MAF requirement 
								last_loc = current_loc;
								int updated_num_var_this_chr=num_geno_variant_array.get(current_chr_index)+1;
								num_geno_variant_array.set(current_chr_index, updated_num_var_this_chr);
								geno_variant_locs_array.get(current_chr_index).add(current_loc);
								geno_G_array.get(current_chr_index).add(alleles);
							}
						}
						line=br.readLine();
					}
				}
				this.num_chr=chr_names_array.size();
				this.chr_names=chr_names_array.toArray(new String[this.num_chr]);
				this.num_geno_variant_M=new int[this.num_chr];
				this.geno_variant_locs=new int[this.num_chr][];
				this.geno_G=new double[this.num_chr][][];
				for(int c=0;c<this.num_chr;c++) {
					this.num_geno_variant_M[c]=num_geno_variant_array.get(c);
					this.geno_variant_locs[c]=new int[this.num_geno_variant_M[c]];
					this.geno_G[c]=new double[this.num_geno_variant_M[c]][];
					for(int m=0;m<this.num_geno_variant_M[c];m++) {
						this.geno_variant_locs[c][m]=geno_variant_locs_array.get(c).get(m);
						this.geno_G[c][m]=geno_G_array.get(c).get(m);
					}	
				}
			}
			else if(this.genotype_file_format.equals("PLINK")) {
				// to be done in the future 
			}else if(this.genotype_file_format.equals("VCF")) {
				// to be done in the future 
			}else {
				System.out.println("Error: Genotype File Format "+
						this.genotype_file_format+" is not supported");
			}
			br.close();
		}catch(Exception e) {e.printStackTrace();}
		
	}
	
	/*
	 * Assuming an 0,1,2 coding of genotype, calculate minor allele frequency. 
	 */
	public static double maf(double[] alleles) {
		double sum = 0;
		for(int n=0;n<alleles.length;n++) {
			sum=sum+alleles[n];
		}
		double maf = sum/(2.0 * alleles.length);
		if(maf>0.5) maf= 1.0 - maf;
		return maf;
	}
	
	/*
	 * set up public HashMap<String, double[][]> gene2filtered_cis_var after readin_genotype_with_maf_ld_filters()
	 */
	public void set_gene2cis_var_map() {
		this.gene2filtered_cis_var = new HashMap<String, double[][]>();
		for(int k=0;k<this.num_gene_K;k++) {
			int[] gene_chr_locs=this.gene_locs.get(gene_names[k]);  //[0]=chr_index; [1]=start; [2]=end.
			int chr_index=gene_chr_locs[0];
			int start=gene_chr_locs[1] - this.cis_var_flanking;
			int end = gene_chr_locs[2] + this.cis_var_flanking;
			// Note that it is OK if start < 0 or end > the last location as the insertion point will ensure this. 
			end = end + this.cis_var_flanking; 
			int start_index_in_geno_G = Arrays.binarySearch(this.geno_variant_locs[chr_index], start);
			if(start_index_in_geno_G<0) { //Arrays.binarySearch returns (-(insertion point) - 1) if not found
				start_index_in_geno_G= -(start_index_in_geno_G + 1);
			}
			int end_index_in_geno_G = Arrays.binarySearch(this.geno_variant_locs[chr_index], end);
			if(end_index_in_geno_G<0) { //Arrays.binarySearch returns (-(insertion point) - 1) if not found
				end_index_in_geno_G= -(start_index_in_geno_G + 1) - 1; 
			}
			double[][] cis_variants_in_gene=new double[end_index_in_geno_G-start_index_in_geno_G+1][];
			for(int m=start_index_in_geno_G; m<=end_index_in_geno_G; m++) { // Note that we don't clone the genotype array. So it doesn't cost memory.
				cis_variants_in_gene[m-start_index_in_geno_G]=this.geno_G[chr_index][m];
			}
			this.gene2filtered_cis_var.put(gene_names[k], cis_variants_in_gene);
		}
	}
	
	/*
	 * Read in Pathway membership from a file
	 * 
	 * Pathway_ID"\t"Gene1,Gene2,...Genen
	 * 
	 * the first column is "\t" separated from the second; and the genes in the second column uses ",". 
	 *
	 * It will generate:
	 *  	public String[] pathway_names;
			public HashMap<String, String[]> pathways
	 *  	public HashMap<String, ArrayList<String>> gene2pathways 
	 */
	public void read_in_pathways() {
		try {
			BufferedReader br= new BufferedReader(new FileReader(this.pathway_file));
			String line=br.readLine();
			while(line!=null) {
				// skip headers
				if(line.startsWith("#"))line=br.readLine(); 
				// parse parameters from the file
				String[] para=line.split("\t"); //para[0]=Pathway_ID; para[1]=list of gene_IDs
				String pathway=para[0];
				String[] genes= para[1].split(",");
				this.pathways2genes.put(pathway, genes); 
				for(int g=0;g<genes.length;g++) {
					if(this.gene2pathways.containsKey(genes[g])) {
						ArrayList<String> pathways=this.gene2pathways.get(genes[g]);
						pathways.add(pathway);
						this.gene2pathways.put(genes[g], pathways);
					}else {
						ArrayList<String> pathways=new ArrayList<String>();
						pathways.add(pathway);
						this.gene2pathways.put(genes[g], pathways);
					}
				}
				line=br.readLine();
			}br.close();			
			this.pathway_names=this.pathways2genes.keySet().toArray(new String[this.pathways2genes.size()]);
		}catch(Exception e) {e.printStackTrace();}
	}
	
	/*
	 * Read in genes_file
	 * the format looks like: gene_ID	chr_num	start_loc	end_loc
	 * columns separated by "\t"
	 * 
	 * It will initiate the following attributes:
			public int num_gene_K;		
			public String[] gene_names;	
			public HashMap<String, Integer> gene_names2index
			public HashMap<String, int[]> gene_locs 
			public boolean[] gene_exp_finalized; 
	 * 
	 * Note that the chr, start, and end are from 1, consistent to the typical GFF files. although in the code, the indexes start at 0. 
	 */
	public void read_in_genes_locs() {
		ArrayList<String> gene_names_list=new ArrayList<String>();
		try {
			BufferedReader br= new BufferedReader(new FileReader(this.genes_file));
			String line=br.readLine();
			while(line!=null) {
				// skip headers
				if(line.startsWith("#"))line=br.readLine(); 
				// parse parameters from the file
				String[] para=line.split("\t"); //para[0]=Gene_ID; para[1]=chr_value; para[2]=start_loc; para[3]=end_loc.
				int[] gene_locations=new int[3];
				gene_locations[0]=this.chr2index.get(para[1]);
				gene_locations[1]=Integer.parseInt(para[2]);
				gene_locations[2]=Integer.parseInt(para[3]);
				gene_names_list.add(para[0]);
				this.gene_locs.put(para[0], gene_locations);
				line=br.readLine();
			}br.close();			
		}catch(Exception e) {e.printStackTrace();}
		this.num_gene_K=gene_names_list.size();
		this.gene_exp_finalized=new boolean[this.num_gene_K];
		Arrays.fill(this.gene_exp_finalized, false); 
		this.gene_names=new String[this.num_gene_K];
		for(int k=0;k<this.num_gene_K;k++) {
			this.gene_names[k]=gene_names_list.get(k);
			this.gene_names2index.put(this.gene_names[k], k);
		}
	}
	
	/* 
	 * set up contributing terms and interacting model of traits and expressions. 
	 * It will initiate 
	 * 		public String[] trait_names;// in the format of T_ID;
			public HashMap<String, Integer> trait_names2index = new HashMap<String, Integer>();
			public boolean[] trait_finalized;
		and generate the lines of terms' contributors and write to:
			this.arch_detailed_file 
		as well as the causality graph corresponding to the arch_detailed_file: 	
			public HashMap<String, ArrayList<String>> causality_graph;
	 */
	public void set_terms_model() {
		this.trait_names=new String[this.num_trait_T]; 
		this.trait_finalized=new boolean[this.num_trait_T];
		Arrays.fill(this.trait_finalized, false);
		for(int t=0;t<this.num_trait_T;t++) {
			this.trait_names[t]=(generator.nextDouble()<=this.trait_binary_proportion?"Binary":"Quantitative")+"_Trait_"+String.format("%04d",t);
			this.trait_names2index.put(this.trait_names[t], t);
		}
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(this.arch_detailed_file));
			// write parameters used to generate this file.
			bw.write("##genotype_file="+this.genotype_file+"\n");
			bw.write("##genes_file="+this.genes_file+"\n");
			bw.write("##pathway_file="+this.pathway_file+"\n");
			bw.write("##max_maf="+this.max_maf+"\n");
			bw.write("##min_maf="+this.min_maf+"\n");
			bw.write("##location_distance="+this.location_distance+"\n");
			bw.write("##gene_contributors={");
			for(int i=0;i<this.gene_contributors.length;i++)
				bw.write(this.gene_contributors[i]+((i==this.gene_contributors.length-1)?"}\n":";"));
			bw.write("##traits_contributors={");
			for(int i=0;i<this.traits_contributors.length;i++)
				bw.write(this.traits_contributors[i]+((i==this.traits_contributors.length-1)?"}\n":","));
			bw.write("##weights_relative_range="+this.weights_relative_range+"\n");
			bw.write("##var_comp_inf_min="+this.traits_var_comp_inf_min+"\n");
			bw.write("##var_comp_inf_max="+this.traits_var_comp_inf_max+"\n");
			bw.write("##var_comp_noise_min="+this.var_comp_noise_min+"\n");
			bw.write("##var_comp_noise_max="+this.var_comp_noise_max+"\n");
			
			bw.write("##cis_variant_numbers_mean="+this.cis_variant_numbers_mean+"\n");
			bw.write("##cis_variant_numbers_range="+this.cis_variant_numbers_range+"\n");
			bw.write("##trans_variant_numbers_mean="+this.trans_variant_numbers_mean+"\n");
			bw.write("##trans_variant_numbers_range="+this.trans_variant_numbers_range+"\n");
			bw.write("##num_contributing_genes_mean="+this.num_contributing_genes_mean+"\n");
			bw.write("##num_contributing_genes_range="+this.num_contributing_genes_range+"\n");
			bw.write("##num_contributing_traits_mean="+this.num_contributing_traits_mean+"\n");
			bw.write("##num_contributing_traits_range="+this.num_contributing_traits_range+"\n");
			bw.write("##gene_models={");
			for(int i=0;i<this.gene_models.length;i++)
				bw.write(this.gene_models[i]+((i==this.gene_models.length-1)?"}\n":";"));
			bw.write("##trait_models={");
			for(int i=0;i<this.trait_models.length;i++)
				bw.write(this.trait_models[i]+((i==this.trait_models.length-1)?"}\n":";"));
			bw.write("##trait_binary_proportion="+this.trait_binary_proportion+"\n");
			// write the column header
			bw.write("#term_ID\tnum_cis\tnum_trans\ttrans_genes\ttrans_exp\ttraits\tWeights\tModel"
					+ "\tInfinitesimal\tNoiseVarianceComponent\n");
			// WRITE GENES FIRST. assign the terms randomly using the probabilistic parameters
			for(int k=0;k<this.num_gene_K;k++) {
				ArrayList<String> causal_terms_of_this_gene = new ArrayList<String>();
				double[] actual_g_weights = new double[this.gene_contri_weights.length]; // [0]=cis;[1]=trans;[2]=other_gene_exp;[3]=trait. 
				//term_ID
				bw.write(this.gene_names[k]+"\t"); 
				// number of cis genetic variants
				if(this.generator.nextDouble()<this.gene_contributors[0]) {  // this gene indeed has [0]=cis contributor(s)
					int cis_number=(int)(cis_variant_numbers_mean+(2*this.generator.nextDouble()-1.0)*(this.cis_variant_numbers_range));
					bw.write((cis_number>=1?cis_number:1)+"\t"); // num_cis written and assign the weight below
					actual_g_weights[0]=this.gene_contri_weights[0]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range*gene_contri_weights[0]);
				}else {
					bw.write("-1\t"); // num_cis is NaN
					actual_g_weights[0]=Double.NaN;
				}
				// number of trans genetic variants
				if(this.generator.nextDouble()<this.gene_contributors[1]) {  // this gene indeed has [1]=trans contributor(s)
					int trans_var_num=(int)(trans_variant_numbers_mean+(2*this.generator.nextDouble()-1.0)*(this.trans_variant_numbers_range));
					bw.write((trans_var_num>=1?trans_var_num:1)+"\t"); // num_trans written and assign the weight below
					actual_g_weights[1]=this.gene_contri_weights[1]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range*gene_contri_weights[1]);
				// genes that the above trans genetic variants will locate. Priority is given to the genes in the same pathway of the focal gene
					int trans_gene_num=(int)(num_contributing_genes_mean+(2*this.generator.nextDouble()-1.0)*(this.num_contributing_genes_range));
					if(trans_gene_num<=1) trans_gene_num=1;
					String[] genes_in_pathway=sample_genes_in_pathways(gene_names[k], trans_gene_num);
					for(int kp=0;kp<genes_in_pathway.length;kp++) {
						causal_terms_of_this_gene.add(genes_in_pathway[kp]+"_Genet");
						bw.write(genes_in_pathway[kp]+((kp==genes_in_pathway.length-1)?"\t":",")); // trans_genes written
					}						
				}else {
					bw.write("-1\tNA\t"); // num_trans is NA, so does trans_genes
					actual_g_weights[1]=Double.NaN;
				}
				// trans gene expressions
				if(this.generator.nextDouble()<this.gene_contributors[2]) {  // this gene indeed has [2]=other gene expression contributor(s)
					int exp_gene_num=(int)(num_contributing_genes_mean+(2*this.generator.nextDouble()-1.0)*(this.num_contributing_genes_range));
					if(exp_gene_num<=1) exp_gene_num=1;
					String[] genes_in_pathway=sample_genes_in_pathways(gene_names[k], exp_gene_num);
					for(int kp=0;kp<genes_in_pathway.length;kp++) {
						causal_terms_of_this_gene.add(genes_in_pathway[kp]+"_Exp");
						bw.write(genes_in_pathway[kp]+((kp==genes_in_pathway.length-1)?"\t":",")); // trans expressions written
					}
					actual_g_weights[2]=this.gene_contri_weights[2]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range*gene_contri_weights[2]);
				}else {
					bw.write("NA\t"); // trans_exp is NA
					actual_g_weights[2]=Double.NaN;
				}
				// traits
				if(this.generator.nextDouble()<this.gene_contributors[3]) {  // this gene indeed has [3]=traits contributor(s)
					int trait_num=(int)(num_contributing_traits_mean+(2*this.generator.nextDouble()-1.0)*(this.num_contributing_traits_range));
					if(trait_num<=1) trait_num=1;
					String[] related_traits=sample_traits(trait_num);
					for(int tp=0;tp<related_traits.length;tp++) {
						causal_terms_of_this_gene.add(related_traits[tp]);
						bw.write(related_traits[tp]+((tp==related_traits.length-1)?"\t":",")); // traits written
					}
					actual_g_weights[3]=this.gene_contri_weights[3]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range*gene_contri_weights[3]);
				}else {
					bw.write("NA\t"); // trait is NA,
					actual_g_weights[3]=Double.NaN;
				}	
				// Weights 
				for(int w=0;w<actual_g_weights.length;w++)
					bw.write(actual_g_weights[w]+((w==actual_g_weights.length-1)?"\t":",")); 
				// Model
				bw.write(sample_model(this.gene_models)+"\t");
				// Infinitesimal
				bw.write("NA\t"); // No Infinitesimal for a gene expression (this column is for traits only).
				// NoiseVarianceComponent
				bw.write(this.generator.nextDouble()*(this.var_comp_noise_max-this.var_comp_noise_min)+this.var_comp_noise_min+"\n");
				// Update Causality Graph
				this.causality_graph.put(gene_names[k], causal_terms_of_this_gene);
			}
			// WRITE TRAITS NEXT
			for(int t=0;t<this.num_trait_T;t++) {
				//term ID
				bw.write(this.trait_names[t]+"\t");
				ArrayList<String> causal_terms_of_this_trait = new ArrayList<String>();
				double[] actual_t_weights = new double[this.traits_contri_weights.length]; // [0]=genetics;[1]=gene_exp;[2]=trait. 
				// there is no cis genetic variants for a term!
				bw.write("-1\t"); // num_cis is NA
				// number of genetic variants
				if(this.generator.nextDouble()<this.traits_contributors[0]) {  // this trait indeed has [0]=genetic contributor(s)
					int trans_var_num=(int)(trans_variant_numbers_mean+(2*this.generator.nextDouble()-1.0)*(this.trans_variant_numbers_range));
					bw.write((trans_var_num>=1?trans_var_num:1)+"\t"); // num_trans written and assign the weight below
					actual_t_weights[0]=this.traits_contri_weights[0]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range * this.traits_contri_weights[0]);
				// genes that the above genetic variants will locate. No consideration of pathways
					int wg_gene_num=(int)(num_contributing_genes_mean+(2*this.generator.nextDouble()-1.0)*(this.num_contributing_genes_range));
					if(wg_gene_num<=1) wg_gene_num=1;
					String[] wg_genes_selected=sample_genes(wg_gene_num);
					for(int kp=0;kp<wg_genes_selected.length;kp++) {
						causal_terms_of_this_trait.add(wg_genes_selected[kp]+"_Genet");
						bw.write(wg_genes_selected[kp]+((kp==wg_genes_selected.length-1)?"\t":",")); // trans_genes written
					}						
				}else {
					bw.write("-1\tNA\t"); // num_trans is NA, so does trans_genes
					actual_t_weights[0]=Double.NaN;
				}
				// gene expressions
				if(this.generator.nextDouble()<this.traits_contributors[1]) {  // this trait indeed has [1]=gene expression contributor(s)
					int exp_gene_num=(int)(num_contributing_genes_mean+(2*this.generator.nextDouble()-1.0)*(this.num_contributing_genes_range));
					if(exp_gene_num<=1) exp_gene_num=1;
					String[] wg_genes_selected=sample_genes(exp_gene_num);
					for(int kp=0;kp<wg_genes_selected.length;kp++) {
						causal_terms_of_this_trait.add(wg_genes_selected[kp]+"+Exp");
						bw.write(wg_genes_selected[kp]+((kp==wg_genes_selected.length-1)?"\t":",")); // gene expressions written
					}
					actual_t_weights[1]=this.traits_contri_weights[1]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range*traits_contri_weights[1]);
				}else {
					bw.write("NA\t"); // trans_exp is NA
					actual_t_weights[1]=Double.NaN;
				}
				// other traits
				if(this.generator.nextDouble()<this.traits_contributors[2]) {  // this trait indeed has [2]=other traits contributor(s)
					int trait_num=(int)(num_contributing_traits_mean+(2*this.generator.nextDouble()-1.0)*(this.num_contributing_traits_range));
					if(trait_num<=1) trait_num=1;
					String[] related_traits=sample_traits(trait_num);
					for(int tp=0;tp<related_traits.length;tp++) {
						causal_terms_of_this_trait.add(related_traits[tp]);
						bw.write(related_traits[tp]+((tp==related_traits.length-1)?"\t":",")); // traits written
					}
					actual_t_weights[2]=this.traits_contri_weights[2]+(2*this.generator.nextDouble()-1.0)*
							(this.weights_relative_range*traits_contri_weights[2]);
				}else {
					bw.write("NA\t"); // trait is NA,
					actual_t_weights[2]=Double.NaN;
				}	
				// Weights 
				for(int w=0;w<actual_t_weights.length;w++)
					bw.write(actual_t_weights[w]+((w==actual_t_weights.length-1)?"\t":",")); 
				// Model
				bw.write(sample_model(this.trait_models)+"\t");
				// Infinitesimal
				if(this.generator.nextDouble()<this.traits_contributors[3]) {  // this trait indeed has [3]=infinitesimal contributor(s)
					bw.write(this.generator.nextDouble()*(this.traits_var_comp_inf_max-this.traits_var_comp_inf_min)+this.traits_var_comp_inf_min+"\t"); 
				}else {
					bw.write("NaN\t"); // Infinitesimal is NaN,
				}
				// NoiseVarianceComponent
				bw.write(this.generator.nextDouble()*(this.var_comp_noise_max-this.var_comp_noise_min)+this.var_comp_noise_min+"\n");
				// Update Causality Graph:
				this.causality_graph.put(this.trait_names[t], causal_terms_of_this_trait);
			}
			bw.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	/*
	 * Sample genes in the pathways related to the focal gene. 
	 * If the focal genes isn't in any pathways, then use the first few pathways in the this.pathway_names
	 */
	public String[] sample_genes_in_pathways(String the_gene_name, int num_genes_needed){
		if(num_genes_needed==0) return null;
		double candidate_ratio = 2.0;  // retain at least twice of genes than the required num_genes_needed.
		double min_candidates=candidate_ratio*num_genes_needed;
		if(this.gene2pathways.containsKey(the_gene_name)) { // this gene is in pathway(s)
			ArrayList<String> pathways = this.gene2pathways.get(the_gene_name);	
			ArrayList<String> genes_in_pathways = new ArrayList<String>();
			for(int i=0;i<pathways.size();i++) { // take the first few pathways until the total number of genes are higher than twice of the num_genes_needed
				if(genes_in_pathways.size()<min_candidates) {
					String[] the_genes=this.pathways2genes.get(pathways.get(i));
					for(int k=0;k<the_genes.length;k++) genes_in_pathways.add(the_genes[k]);
				}else break;
			}
			if(genes_in_pathways.size()<=num_genes_needed) { // no enough candidates, just return all of them
				return genes_in_pathways.toArray(new String[genes_in_pathways.size()]);
			}else { // more candidates than num_genes_needed. So select randomly.
				String[] selected_genes=new String[num_genes_needed];
				for(int k=0;k<num_genes_needed;k++) {
					int index=this.generator.nextInt(genes_in_pathways.size());
					selected_genes[k]=genes_in_pathways.get(index);
					genes_in_pathways.remove(index);  // remove the selected.
				}
				return selected_genes;
			}
		}else {  // this gene is not in any pathways: sample a list of genes from randomly selected pathways 
			ArrayList<String> genes_in_pathways = new ArrayList<String>();
			for(int i=0;i<this.pathway_names.length;i++) { // take the first few pathways until the total number of genes are higher than twice of the num_genes_needed
				if(genes_in_pathways.size()<min_candidates) {
					String[] the_genes=this.pathways2genes.get(this.pathway_names[i]);
					for(int k=0;k<the_genes.length;k++) genes_in_pathways.add(the_genes[k]);
				}else break;
			}
			if(genes_in_pathways.size()<=num_genes_needed) { // no enough candidates, just return all of them
				return genes_in_pathways.toArray(new String[genes_in_pathways.size()]);
			}else { // more candidates than num_genes_needed. So select randomly.
				String[] selected_genes=new String[num_genes_needed];
				for(int k=0;k<num_genes_needed;k++) {
					int index=this.generator.nextInt(genes_in_pathways.size());
					selected_genes[k]=genes_in_pathways.get(index);
					genes_in_pathways.remove(index);  // remove the selected.
				}
				return selected_genes;
			}
		}
	}
	
	/*
	 * Randomly select several traits from this.trait_names[].
	 */
	public String[] sample_traits(int number_traits_needed) {
		String[] selected=new String[number_traits_needed];
		for(int t=0;t<number_traits_needed;t++) {
			selected[t]=this.trait_names[this.generator.nextInt(this.num_trait_T)];
			// could select the same trait twice, which is OK.
		}
		return selected;
	}
	
	/*
	 * Randomly select several genes randomly without considering pathways from this.trait_names[].
	 */
	public String[] sample_genes(int number_genes_needed) {
		String[] selected=new String[number_genes_needed];
		for(int k=0;k<number_genes_needed;k++) {
			selected[k]=this.gene_names[this.generator.nextInt(this.num_gene_K)];
			// could select the same gene twice, which is OK.
		}
		return selected;
	}
	
	/*
	 * Randomly select from SpecificModels.supported_models={"additive", "epistatic", "compensatory", "heterogenous", "compound"} 
	 * based on their probabilities in frequency[].
	 */
	public String sample_model(double[] frequency) {
		double the_random = this.generator.nextDouble();
		double sum=0;
		for(int i=0;i<frequency.length;i++) sum+=frequency[i];
		if(Double.compare(sum, 1.0)!=0) {
			System.out.println("Error: Sum of the frequency array not 1.0 in the function String sample_model(double[] frequency)");
			return null;
		}
		double current_cut=0.0;
		for(int i=0;i<frequency.length;i++) {
			if(the_random>=current_cut && the_random<current_cut+frequency[i])
				return SpecificModels.supported_models[i];
			current_cut+=frequency[i];	
		}
		return null; // this will not be executed as one of the return in the for loop must be executed. 
	}
	
	/*
	 * Randomly select from a chromsome.
	 * The chance of landing a chromosome is proportional to its number of genetic variants 
	 * Return the index of the selected chromosome 
	 */
	public int sample_chr() {
		double the_random = this.generator.nextDouble();
		double[] frequency = SpecificModels.to_prob_distribution(this.num_geno_variant_M.clone());
		double sum=0;
		for(int i=0;i<this.num_chr;i++) sum+=frequency[i];
		if(Double.compare(sum, 1.0)!=0) {
			System.out.println("Error: Sum of the frequency array not 1.0 in the function String sample_chr()");
			return -1;
		}
		double current_cut=0.0;
		for(int i=0;i<frequency.length;i++) {
			if(the_random>=current_cut && the_random<current_cut+frequency[i])
				return i;
			current_cut+=frequency[i];	
		}
		return -1; // this will not be executed as one of the return in the for loop must be executed. 
	}
	
	/* 	
	 * Calculate 
	 * 	public double[][] corr_KxK;	// correlations between gene expressions
	 * 	public double[][] corr_TxT;	// correlations between traits 
	 *	public double[][] corr_TxK;	// correlations between traits and expressions
		Note that this function may be used in any time, however its current use in OmeSim is 
		BEFORE adding the noise and infinitesmal terms to reflect genuine correlations 
	 */
	
	public void calculate_correlations() { 
		// assign public double[][] corr_KxK;	// correlations between gene expressions
		this.corr_KxK=new double[this.num_gene_K][]; // a triangle matrix 
		for(int k1=0; k1<this.num_gene_K; k1++) {
			this.corr_KxK[k1]=new double[this.num_gene_K-k1];
			for(int k2=k1; k2<this.num_gene_K; k2++) {
				this.corr_KxK[k1][k2-k1]=SpecificModels.correlation(this.exp_Z[k1], this.exp_Z[k2]);
			}
		}		
		// assign public double[][] corr_TxT;	// correlations between traits
		this.corr_TxT=new double[this.num_trait_T][]; // a triangle matrix 
		for(int t1=0; t1<this.num_trait_T; t1++) {
			this.corr_TxT[t1]=new double[this.num_trait_T-t1];
			for(int t2=t1; t2<this.num_trait_T; t2++) {
				this.corr_TxT[t1][t2-t1]=SpecificModels.correlation(this.trait_Y[t1], this.trait_Y[t2]);
			}
		}	
		// assign public double[][] corr_TxK;	// correlations between traits and expressions
		this.corr_TxK=new double[this.num_gene_K][this.num_trait_T];
		for(int k1=0; k1<this.num_gene_K; k1++) {
			for(int t1=0; t1<this.num_trait_T; t1++) {
				this.corr_TxK[t1][k1]=SpecificModels.correlation(this.trait_Y[t1], this.exp_Z[k1]);
			}
		}
	}
	
	/*
	 * Calculate 
	 *  public double[][] asso_MxK;	// associations between individual genetic variants and expressions 
	 * 	public double[][] asso_MxT;	// associations between individual genetic variants and traits 
	 	Note that this function may be used in any time, however its current use in OmeSim is 
		BEFORE adding the noise and infinitesimal terms to reflect genuine correlations 
	 */
	public void calculate_associations(CausalTerms[] full_terms) { //TODO
		//this.asso_MxK
	}
	/*
	 * Output generated data to files! TODO
	 * the following data will be written to files:
	 */
	public void output_files(String output_file_folder) {
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(""));
			bw.write("#gene_contri_weights={");
			for(int i=0;i<this.gene_contri_weights.length;i++)
				bw.write(this.gene_contri_weights[i]+((i==this.gene_contri_weights.length-1)?"}\n":","));
			bw.write("#traits_contri_weights={");
			for(int i=0;i<this.traits_contri_weights.length;i++)
				bw.write(this.traits_contri_weights[i]+((i==this.traits_contri_weights.length-1)?"}\n":";"));
			bw.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	
	
}
