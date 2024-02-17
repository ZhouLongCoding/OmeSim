package omesim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

/*
 * Preparing files and call external scripts to generate plots for
 * 1. Expression / Traits values
 * 2. Correlation files 
 * 3. Causality graphs
 * 
 * Both overall and pathways-based figures will be generated.
 * 
 */
public class Visualization {
	public MainFrame main_frame;
	public String input_data_folder;
	public String figure_output_folder;
	public String corr_R_code_file;
	public String data_R_code_file;
	public String causality_R_code_file;
	public String Rscript_binary_path;
	public int trait_causality_degree;
	public double edge_num_ratio;
//	public boolean large_figures_needed;
	
	/*
	 * Constructor setting the fields and making the figures 
	 */
	public Visualization(MainFrame main_frame) {
		this.main_frame=main_frame;
		this.input_data_folder=main_frame.output_folder;
		this.figure_output_folder=main_frame.output_folder+"/figures/";
		File figures_folder=new File(this.figure_output_folder);
		if(!figures_folder.exists()) {
			figures_folder.mkdir();
		}
		this.corr_R_code_file = main_frame.visualization_code_folder+"/corr_heatmap.R";
		this.data_R_code_file= main_frame.visualization_code_folder+"/data_heatmap.R";
		this.causality_R_code_file= main_frame.visualization_code_folder+"/causality_network.R";
		this.Rscript_binary_path=main_frame.Rscript_binary_path;
		//this.large_figures_needed=main_frame.large_figures_needed;
		this.trait_causality_degree=this.main_frame.trait_causality_degree;
		this.edge_num_ratio=this.main_frame.edge_num_ratio;
	}
	
	/*
	 * make all figures (including overall and pathway-based.
	 */
	public void make_all_figures(boolean make_full) {
		System.out.println("Started Drawing graphics.");
		if(make_full){
			// make 5 plots for full figures first
			this.make_full_pic_plots();
		}
		
		// split data into pathways or traits
		this.generate_pathway_files_for_data_and_correlations();
		int[] pathway_edge_nums = this.make_causality_pathway_files();
		int[][] trait_edge_nums = this.make_causality_trait_files(this.trait_causality_degree);
		
		// draw figures for each pathway or trait
		this.draw_pathway_figures_for_data_and_correlations();
		this.draw_causality_figures(pathway_edge_nums, trait_edge_nums, this.trait_causality_degree);
		System.out.println("Graphics Done.");
	}
	
	/*
	 * Split large files into smaller ones based on pathways. Three files in each pathway:
	 * 
	 * Gene expression data
	 * Correlation: Exp x Exp
	 * Correlation: Exp x Trait  
	 */
	public void generate_pathway_files_for_data_and_correlations() {
		try {			
			// split related files to pathways. 
			// iteration on each pathways
			for(int p=0;p<this.main_frame.pathway_names.length;p++) {
				String folder_pathway_string=this.figure_output_folder+"/"+this.main_frame.pathway_names[p];
				File folder_pathway=new File(folder_pathway_string);
				if(!folder_pathway.exists()) folder_pathway.mkdir();
				String[] genes_pathw= this.main_frame.pathways2genes.get(this.main_frame.pathway_names[p]);
				// extract data 
				int[] gene_indexes_full=new int[genes_pathw.length];
				for(int k_pathw=0;k_pathw<genes_pathw.length;k_pathw++) { // note that g_pathw is the index in genes, not the full datasets. 
					// get the index in the main datasets for each gene. 
					gene_indexes_full[k_pathw]=this.main_frame.gene_names2index.get(genes_pathw[k_pathw]);
				}
				
				// Get gene expression data and write to file
				String exp_pathw_csv_file=folder_pathway+"/"+
						this.main_frame.pathway_names[p]+".expressions.csv";
				BufferedWriter bw= new BufferedWriter(new FileWriter(exp_pathw_csv_file));
				bw.write("#Gene_ID");
				for(int i=0;i<this.main_frame.num_subj_N;i++) {
					bw.write(","+this.main_frame.subj_names[i]);
				}
				bw.write("\n"); // finished writing the header line, which is the sample names
				
				for(int k_pathw=0;k_pathw<genes_pathw.length;k_pathw++) {
					bw.write(genes_pathw[k_pathw]); // expressions of one gene
					for(int i=0;i<this.main_frame.num_subj_N;i++) {
						bw.write(","+String.format("%.4f",this.main_frame.exp_Z[gene_indexes_full[k_pathw]][i]));
					}
					bw.write("\n");
				}
				bw.close();
//				// plot expression data
//				String[] R_args_exp= {exp_pathw_csv_file, 
//						folder_pathway_string+"/"+this.main_frame.pathway_names[p]+".expressions.png",
//						"Expressions with Clustering in "+this.main_frame.pathway_names[p]};
//				this.make_a_plot(folder_pathway_string, this.data_R_code_file, R_args_exp);
					
				// Get correlation: Exp x Exp and write to file
				String e_by_e_csv_file=folder_pathway+"/"
						+this.main_frame.pathway_names[p]+".corr_exp_by_exp.csv";
				bw= new BufferedWriter(new FileWriter(e_by_e_csv_file));
				bw.write("#Gene_ID");
				for(int k_pathw=0;k_pathw<genes_pathw.length;k_pathw++) {
					bw.write(","+genes_pathw[k_pathw]);
							//+((k==this.num_gene_K-1)?"\n":",")); 
				}bw.write("\n");
				// finished writing the header line
				// to facilitate the clustering within heatmap, we output a rectangle instead of a triangle 
				for(int k1=0; k1<genes_pathw.length; k1++) {
					bw.write(genes_pathw[k1]); 
					for(int k2=0; k2<genes_pathw.length; k2++) {
							bw.write(","+String.format("%.4f",
									this.main_frame.corr_KxK[gene_indexes_full[k1]][gene_indexes_full[k2]]));
					//				+((k2==this.num_gene_K-1)?"\n":",")); 					
					}bw.write("\n");
				}	
				bw.close();
//				// Plot correlation: Exp x Exp
//				String[] R_args_exp_by_exp= {e_by_e_csv_file, 
//						folder_pathway_string+"/"+this.main_frame.pathway_names[p]+".corr_exp_by_exp.png",
//						"Correlations between Genes in "+this.main_frame.pathway_names[p]};
//				this.make_a_plot(folder_pathway_string, this.corr_R_code_file, R_args_exp_by_exp);
			
				// Get correlation: Trait x Exp and write to file
				String t_by_e_csv_file = folder_pathway+"/"
						+this.main_frame.pathway_names[p]+".trait_by_exp.csv";
				bw= new BufferedWriter(new FileWriter(t_by_e_csv_file));
				
				bw.write("#Gene_ID");
				for(int t=0;t<this.main_frame.num_trait_T;t++) {
					bw.write(","+this.main_frame.trait_names[t]); 
				}
				bw.write("\n");  // finished writing the header line 
				for(int k1=0; k1<genes_pathw.length; k1++) {
					bw.write(genes_pathw[k1]);
					for(int t1=0; t1<this.main_frame.num_trait_T; t1++) {
						bw.write(","+String.format("%.4f",
								this.main_frame.corr_TxK[t1][gene_indexes_full[k1]]));
					}bw.write("\n");
				}				
				bw.close();
				// Plot correlation: Trait x Exp
//				String[] R_args_trait_by_exp= {t_by_e_csv_file, 
//						folder_pathway_string+"/"+this.main_frame.pathway_names[p]+".trait_by_exp.png",
//						"Correlations between Traits and Genes in "+this.main_frame.pathway_names[p]};
//				this.make_a_plot(folder_pathway_string, this.corr_R_code_file, R_args_trait_by_exp);
			}
			 
			//
			System.out.println("Finished generating data/correlation files for individual pathways");
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public void draw_pathway_figures_for_data_and_correlations() {
		try {			
			// split related files to pathways. 
			// iteration on each pathways
			for(int p=0;p<this.main_frame.pathway_names.length;p++) {
				String folder_pathway=this.figure_output_folder+"/"+this.main_frame.pathway_names[p];
				// plot expression data
				String exp_pathw_csv_file=folder_pathway+"/"+
						this.main_frame.pathway_names[p]+".expressions.csv";
				String[] R_args_exp= {exp_pathw_csv_file, 
						folder_pathway+"/"+this.main_frame.pathway_names[p]+".expressions.png",
						"Expressions with Clustering in "+this.main_frame.pathway_names[p]};
				this.make_a_plot(folder_pathway, this.data_R_code_file, R_args_exp);
		
				// Plot correlation: Exp x Exp
				String e_by_e_csv_file=folder_pathway+"/"
						+this.main_frame.pathway_names[p]+".corr_exp_by_exp.csv";
				String[] R_args_exp_by_exp= {e_by_e_csv_file, 
						folder_pathway+"/"+this.main_frame.pathway_names[p]+".corr_exp_by_exp.png",
						"Correlations between Genes in "+this.main_frame.pathway_names[p]};
				this.make_a_plot(folder_pathway, this.corr_R_code_file, R_args_exp_by_exp);
						
				// Plot correlation: Trait x Exp
				String t_by_e_csv_file = folder_pathway+"/"
						+this.main_frame.pathway_names[p]+".trait_by_exp.csv";
				String[] R_args_trait_by_exp= {t_by_e_csv_file, 
						folder_pathway+"/"+this.main_frame.pathway_names[p]+".trait_by_exp.png",
						"Correlations between Traits and Genes in "+this.main_frame.pathway_names[p]};
				this.make_a_plot(folder_pathway, this.corr_R_code_file, R_args_trait_by_exp);
			}
			System.out.println("Finished drawing figures of data and correlations for individual pathways");
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/*
	 * generate pathway files for causality plots using 'igraph' in R
	 * 
	 * For each pathway, a causality network will be plotted
	 * 		Note that all edges having the members of a pathway as a target will be extracted to plot the 
	 * 		sub-causality-graph of that pathway.  
	 * 
	 *  
	 */
	public int[] make_causality_pathway_files() {
		// count the number of edges in this graph to set up the size of the figure
		int[] pathway_edge_num=new int[this.main_frame.pathway_names.length]; 
		try {
			// create files for pathways
			for(int p=0;p<this.main_frame.pathway_names.length;p++) {
				// create the folder and the file
				String folder_pathway_string=this.figure_output_folder+"/"+this.main_frame.pathway_names[p];
				File folder_pathway=new File(folder_pathway_string);
				if(!folder_pathway.exists()) folder_pathway.mkdir();
				BufferedWriter bw_pathway=new BufferedWriter(new FileWriter(folder_pathway+"/"+
						this.main_frame.pathway_names[p]+".network.csv"));
				bw_pathway.write("#Source,Target\n");
				// all genes in this pathway
				String[] genes_in_pathw=this.main_frame.pathways2genes.get(this.main_frame.pathway_names[p]);
				// for each gene, process all sources for this gene.
				
				for(int g=0;g<genes_in_pathw.length;g++) {
					ArrayList<String> sources=this.main_frame.causality_graph.get(genes_in_pathw[g]);
					for(int source_index=0;source_index<sources.size();source_index++) {
						bw_pathway.write(sources.get(source_index)+","+genes_in_pathw[g]+"\n");
						pathway_edge_num[p]++;
					}
				}
				bw_pathway.close();
			}			
		}catch (Exception e){
			e.printStackTrace();
		}
		return pathway_edge_num;
	}
	
	/*
	 * generate trait files for causality plots using 'igraph' in R
	 * 
	 *  For each trait, several causality networks based on different degrees will be plotted.
	 * 		Here the "degree" means the number of edges through which a node is linked to the trait.
	 */
	public int[][] make_causality_trait_files(int max_network_degree) {
		int[][] trait_edge_nums=new int[this.main_frame.num_trait_T][max_network_degree]; 
		try {
	// create files for traits
			for(int t=0;t<this.main_frame.num_trait_T;t++) {
				// create the folder and the file
				String folder_trait_string=this.figure_output_folder+"/"+this.main_frame.trait_names[t];
				File folder_trait=new File(folder_trait_string);
				if(!folder_trait.exists()) folder_trait.mkdir();
				BufferedWriter[] bw_trait=new BufferedWriter[max_network_degree];
				// each trait has multiple networks for each d\in{1,2,...,max_network_degree}
				for(int d=0;d<max_network_degree;d++) {
					bw_trait[d]=new BufferedWriter(new FileWriter(folder_trait+"/"+
						this.main_frame.trait_names[t]+".d"+(d+1)+".csv"));
					bw_trait[d].write("#Source,Target\n");
				}
				// write the trait-related edges into corresponding files
				ArrayList<String> target_set=new ArrayList<String>();
				ArrayList<String> source_set=new ArrayList<String>();
				target_set.add(this.main_frame.trait_names[t]);
				for(int d=0;d<max_network_degree;d++) {
					for(int target_index=0;target_index<target_set.size();target_index++) {
						ArrayList<String> sources=this.main_frame.causality_graph.get(target_set.get(target_index));
						for(int so_index=0;so_index<sources.size();so_index++) {
							// save this item in the source_set for the next round
							source_set.add(sources.get(so_index));
							// write the files that allow this degree
							for(int d_file=d;d_file<max_network_degree;d_file++) {
								bw_trait[d_file].write(sources.get(so_index)+","+target_set.get(target_index)+"\n");
								trait_edge_nums[t][d_file]++;
							}
						}
					}
					// finished writing this degree. Put the current sources as the targets for the next round
					target_set.clear();
					for(int so_set_index=0;so_set_index<source_set.size();so_set_index++) {
						target_set.add(source_set.get(so_set_index));
					}
					// clear the sources for the next round
					source_set.clear();
				}
				for(int d=0;d<max_network_degree;d++) {
					bw_trait[d].close();
				}
			}	
		}catch (Exception e){
			e.printStackTrace();
		}
		return trait_edge_nums;
	}
	
	public void draw_causality_figures(int[] pathway_edge_nums, int[][] trait_edge_nums, 
			int max_network_degree) {
		try {
			// draw figures for pathways
			for(int p=0;p<this.main_frame.pathway_names.length;p++) {
				String folder_pathway=this.figure_output_folder+"/"+this.main_frame.pathway_names[p];
				// plot the sub-causality graph for this pathway 
				int png_size=(int)(edge_num_ratio*Math.pow(pathway_edge_nums[p],0.4));
				String[] R_args_p_causality= {
						folder_pathway+"/"+this.main_frame.pathway_names[p]+".network.csv", 
						folder_pathway+"/"+this.main_frame.pathway_names[p]+".network.png", 
						png_size+"", 
						"Causality Graph for "+this.main_frame.pathway_names[p]};
				this.make_a_plot(folder_pathway, this.causality_R_code_file, R_args_p_causality);
			}			
			System.out.println("Finished drawing figures of causality graphs for individual pathways");
			// draw figures for traits
			for(int t=0;t<this.main_frame.num_trait_T;t++) {
				String folder_trait=this.figure_output_folder+"/"+this.main_frame.trait_names[t];
				for(int d=0;d<max_network_degree;d++) {
					// plot the sub-causality graph for this pathway and this degree
					int png_size=(int)(edge_num_ratio*Math.pow(trait_edge_nums[t][d],0.4));
					String[] R_args_p_causality= {
							folder_trait+"/"+this.main_frame.trait_names[t]+".d"+(d+1)+".csv",
							folder_trait+"/"+this.main_frame.trait_names[t]+".d"+(d+1)+".png",
							png_size+"",
							"Causality Graph for "+this.main_frame.trait_names[t]+": Degree="+(d+1)};
					this.make_a_plot(folder_trait, this.causality_R_code_file, R_args_p_causality);
				}				
			}	
			System.out.println("Finished drawing figures of causality graphs for individual traits");
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/*
	 * Using existing R script to make one figure. 
	 * 
	 * By defualt 
	 * 	String Rscript_binary_path="/usr/local/bin/Rscript"
	 * 	String R_code_file = the R code for correlation or data heatmap 
	 * 	String R_args[] = {data_file, image_file, image_title}
	 * 
	 */
	public void make_a_plot(String working_dir, String R_code_file, String R_args[]) {		
		//System.out.println("Try command line R");
		try {
		    // add the Rscript path and R code file and its arguments to an array.
			String[] command=new String[R_args.length+2];
			command[0]=this.Rscript_binary_path;
			command[1]=R_code_file;
			for(int i=0;i<R_args.length;i++) command[i+2]=R_args[i];
		    ProcessBuilder processBuilder = new ProcessBuilder(command); 
		    // Setting the directory from which the R code will be executed 
		    processBuilder.directory(new File(working_dir));
		    // Execute the R code 
            Process process = processBuilder.start();
//          // Read the output of the R script
//            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
//            String line;
//            while ((line = reader.readLine()) != null) {
//                System.out.println(line);
//            }
            
            // Check for script execution errors
            int exitCode = process.waitFor();
            if (exitCode == 0) {
                System.out.println(R_code_file+" executed successfully for "+R_args[0]);
            } else {
                System.err.println(R_code_file+" execution failed for "+R_args[0]);
            }
		} catch (IOException | InterruptedException e) {
		    e.printStackTrace();
		}
	}
	
	/*
	 * Make plots for five full datasets, including 
	 * 	2 files for traits and expressions data
	 * 	3 files for correlations 
	 * 
	 * String R_args_XYZ[] = {data_file, image_file, image_title}
	 */
	public void make_full_pic_plots() {
		// expressions.csv 
		String[] R_args_exp= {this.input_data_folder+"/expressions.csv", this.figure_output_folder+"/expressions.png", "Expressions with Clustering"};
		this.make_a_plot(this.figure_output_folder, this.data_R_code_file, R_args_exp);
		// traits.csv
		String[] R_args_traits= {this.input_data_folder+"/traits.csv", this.figure_output_folder+"/traits.png", "Traits with Clustering"};
		this.make_a_plot(this.figure_output_folder, this.data_R_code_file, R_args_traits);
		// corr_exp_by_exp.csv 
		String[] R_args_exp_by_exp= {this.input_data_folder+"/corr_exp_by_exp.csv", this.figure_output_folder+"/corr_exp_by_exp.png", "Correlation between Expressions"};
		this.make_a_plot(this.figure_output_folder, this.corr_R_code_file, R_args_exp_by_exp);
		// corr_trait_by_trait.csv 
		String[] R_args_trait_by_trait= {this.input_data_folder+"/corr_trait_by_trait.csv", this.figure_output_folder+"/corr_trait_by_trait.png", "Correlation between Traits"};
		this.make_a_plot(this.figure_output_folder, this.corr_R_code_file, R_args_trait_by_trait);
		// corr_trait_by_exp.csv 
		String[] R_args_trait_by_exp= {this.input_data_folder+"/corr_trait_by_exp.csv", this.figure_output_folder+"/corr_trait_by_exp.png", "Correlation between Traits and Expressions"};
		this.make_a_plot(this.figure_output_folder, this.corr_R_code_file, R_args_trait_by_exp);
		System.out.println("Finsihed generating full plots.");
	}

//	public static void main(String[] args) {
//		String folder="/Users/quanlong/SH Corp Dropbox/quan long/Qingrun/Writing/Manuscripts/OmeSim/GWAS/";
//		String rScriptPath=folder+"visualization_scripts_folder/corr_heatmap.R";
//				//+ "corr_trait_by_trait.csv corr_trait_by_trait.silly.png TITLE";
//		String[] r_args = {"corr_trait_by_trait.csv", "corr_trait_by_trait.silly.png", "TITLE"};
//
//	}
}
