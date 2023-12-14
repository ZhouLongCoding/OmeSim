package omesim;

import java.util.ArrayList;
import java.util.HashMap;

/*
 * Functions to identify the causal relationship of a triple {gene(DNA), expression, traits} 
 * based on generated causal structure: MainFrame:: HashMap<String, ArrayList<String>> causality_graph
 * 
 * Indirect causality of maximal one extra link is allowed: A => B => C will build up the causality of A => C.
 * However, we do not consider more than one extra link for the moment (2023 July).  
 */
public class ArchitectureChecker {
	
	public MainFrame main_frame;
	public HashMap<String, ArrayList<String>> causality_graph;
	
	public ArchitectureChecker(MainFrame main_frame) {
		this.main_frame= main_frame;
		this.causality_graph=main_frame.causality_graph;
	}
	
	/*
	 * Genetic variations are the cause of transcriptomic variations, 
	 * which in turn causes phenotypic changes. 
	 * Specifically, we have Z ~ f_c(G) + \epsilon, and Y ~ g_c(G, Z) + \epsilon. 
	 * The subindex “c” standards for “causality”. 
	 */
	public boolean causality(String variants_gene_ID, String exp_gene_ID, String trait_ID) {
		variants_gene_ID = variants_gene_ID+"_Genet";
		exp_gene_ID = exp_gene_ID+"_Exp";
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
		variants_gene_ID = variants_gene_ID+"_Genet";
		exp_gene_ID = exp_gene_ID+"_Exp";
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
		variants_gene_ID = variants_gene_ID+"_Genet";
		exp_gene_ID = exp_gene_ID+"_Exp";
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
		for(int m=0;m<potential_mediators.size();m++) {
			ArrayList<String> potential_sources = this.causality_graph.get(potential_mediators.get(m));
			for(int k=0;k<potential_sources.size();k++) {
				if(potential_sources.get(k).equals(source)) return true;
			}			
		}return false;
	}
}
