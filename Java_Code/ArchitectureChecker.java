package omesim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

/*
 * 	Functions to identify the causal relationship of a triple {gene(DNA), expression, traits} 
 *  based on generated causal structure: MainFrame:: HashMap<String, ArrayList<String>> causality_graph
 * 
 * 	Indirect causality of maximal one extra link is allowed: A => B => C will build up the causality of A => C.
 * 	However, we do not consider more than one extra link for the moment (2023 July). 
 * 
 * 	For a triple {variants_gene_ID, exp_gene_ID, trait_ID}, we have the following relationship:
 * 		causality:  variants_gene_ID => exp_gene_ID => trait_ID
 * 		pleiotropy: variants_gene_ID => exp_gene_ID; variants_gene_ID => trait_ID
 * 		reverse_causality: variants_gene_ID => trait_ID => exp_gene_ID 
 * 
 *  Note that if one put another category of term (e.g., put a gene_ID in the third place that is supposed for a 
 *  trait_ID, the program will still work, although the biological interpretation may be different.    
 */
public class ArchitectureChecker {
	
	public MainFrame main_frame;
	public HashMap<String, ArrayList<String>> causality_graph;
	
	public ArchitectureChecker(MainFrame main_frame) {
		this.main_frame= main_frame;
		this.causality_graph=main_frame.causality_graph;
	}
	
	/*
	 * load the causality graph from a file, as written by the 
	 * causality graph part of MainFrame.output_files()
	 * target,source1,source2,...
	 */
	public ArchitectureChecker(String causal_graph_file) {
		this.main_frame=null;
		this.causality_graph=new HashMap<String, ArrayList<String>>();
		try {
			BufferedReader br= new BufferedReader(new FileReader(causal_graph_file));
			String line=br.readLine();
			while(line.startsWith("#"))line=br.readLine();
			while(line!=null) {
				String[] ids=line.split(",");
				String target=ids[0];
				if(ids.length>1) { // there are indeed sources for this target
					ArrayList<String> sources=new ArrayList<String>();
					for(int s=1;s<ids.length;s++)sources.add(ids[s]);
					this.causality_graph.put(target, sources);
				}
				line=br.readLine();
			}br.close(); 
		}catch(Exception e) {e.printStackTrace();}
	}
	/*
	 * Check the causal relationship for one triple.
	 */
	public String causal_check(String variants_gene_ID, String exp_gene_ID, String trait_ID) {
		if(causality(variants_gene_ID, exp_gene_ID, trait_ID)) {
			return "causality: "+variants_gene_ID+" => "+exp_gene_ID+" => "+trait_ID;
		}else if(pleiotropy(variants_gene_ID, exp_gene_ID, trait_ID)) {
			return "pleiotropy: "+variants_gene_ID+" => "+exp_gene_ID+"; "+variants_gene_ID+" => "+trait_ID;
		}else if(reverse_causality(variants_gene_ID, exp_gene_ID, trait_ID)) {
			return "reverse_causality: "+variants_gene_ID+" => "+trait_ID+" => "+exp_gene_ID;
		}else {
			return "No-known-relaion: "+variants_gene_ID+" "+exp_gene_ID+" "+trait_ID;
		}
	}
	
	/*
	 * check the causal relationship for multiple triples in a file, each line of the triple file shuold be:
	 * variants_gene_ID,exp_gene_ID,trait_ID
	 * 
	 * The outcome will be their relationship causal relationship. 
	 */
	public void causal_check_file(String input_triple_lines, String output_relation_file) {
		try {
			BufferedReader br= new BufferedReader(new FileReader(input_triple_lines));
			BufferedWriter bw= new BufferedWriter(new FileWriter(output_relation_file));
			String line=br.readLine();
			while(line.startsWith("#"))line=br.readLine();
			while(line!=null) {
				String[] ids=line.split(",");
				bw.write(causal_check(ids[0], ids[1], ids[2]));
				line=br.readLine();
			}br.close(); bw.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	/*
	 * Genetic variations are the cause of transcriptomic variations, 
	 * which in turn causes phenotypic changes. 
	 * Specifically, we have Z ~ f_c(G) + \epsilon, and Y ~ g_c(G, Z) + \epsilon. 
	 * The subindex “c” standards for “causality”. 
	 */
	public boolean causality(String variants_gene_ID, String exp_gene_ID, String trait_ID) {
		return (direct_cause(variants_gene_ID, exp_gene_ID) || indirect_cause(variants_gene_ID, exp_gene_ID)) &&
				(direct_cause(exp_gene_ID, trait_ID) || indirect_cause(exp_gene_ID, trait_ID));
	}
	
	/*
	 * Genetic variations cause the transcriptomic variations, 
	 * and phenotypic changes via independent models. 
	 * Specifically, we have Z ~ f_p(G) + \epsilon, and Y ~ g_p(G) + \epsilon. 
	 * The subindex “c” standards for “pleiotropy”. Please note that, here “independent” 
	 * means the functions f_p(.) and g_p(.) are independent. Since both Z and Y are 
	 * dependent on G, they are still correlated.
	 */
	public boolean pleiotropy(String variants_gene_ID, String exp_gene_ID, String trait_ID) {
		return (direct_cause(variants_gene_ID, exp_gene_ID) || indirect_cause(variants_gene_ID, exp_gene_ID)) &&
				(direct_cause(variants_gene_ID, trait_ID) || indirect_cause(variants_gene_ID, trait_ID));
	}
	
	/*
	 * Genetic variations are the cause of phenotypic changes, 
	 * which in turn causes transcriptomic variations. 
	 * Specifically, we have Y ~ f_r(G) + \epsilon, and Z ~ g_r(G, Y) + \epsilon. 
	 * The subindex “r” standards for “reverse-causality”. 
	 */
	public boolean reverse_causality(String variants_gene_ID, String exp_gene_ID, String trait_ID) {
		return (direct_cause(variants_gene_ID, trait_ID) || indirect_cause(variants_gene_ID, trait_ID)) &&
				(direct_cause(trait_ID, exp_gene_ID) || indirect_cause(trait_ID, exp_gene_ID));
	}
	
	/*
	 * directly causal. 
	 */
	public boolean direct_cause(String source, String target) {
		ArrayList<String> potential_sources= this.causality_graph.get(target);
		for(int k=0;k<potential_sources.size();k++) {
			if(potential_sources.get(k).equals(source)) return true;
		}return false;
	}
	
	/*
	 * indirectly causal based on a mediator  
	 */
	public boolean indirect_cause(String source, String target) {
		ArrayList<String> potential_mediators= this.causality_graph.get(target);
		for(int medator_index=0;medator_index<potential_mediators.size();medator_index++) {
			ArrayList<String> potential_sources = this.causality_graph.get(potential_mediators.get(medator_index));
			for(int k=0;k<potential_sources.size();k++) {
				if(potential_sources.get(k).equals(source)) return true;
			}			
		}return false;
	}
	
	
}
