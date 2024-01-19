package omesim;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

//import org.apache.commons.math3.distribution.NormalDistribution;

/*
 * For all input data, in this class, denoted as X[][], 
 * the first dimension will be sample size, 
 * and the second is the number of terms 
 * 
 * Note that X[][] could be a mixture of difference omics, if the users prefer. 
 * Then a re-scale by function "standardization(double[] data, double weight)" 
 * may to be carried out for some terms with substantial different scale. 
 * 
 * For logical operations (AND/OR/XOR), <=0 means false. When it is raw genotype, it is 0,1,2.
 * Or if it is expression or rep-learned genotype, its value is centralized by 0,
 * and then <0 means false.  
 * 
 */
public class SpecificModels {
	
	public static final String _additive="additive";
	public static final String _epistatic="epistatic";
	public static final String _compensatory="compensatory";
	public static final String _heterogenous="heterogenous";
	public static final String _compound="compound";
	
	public static final int random_seed =1;
	public static final double default_weight=1.0;
	public static final Random generator = new Random(random_seed);
	
	public static final String[] supported_models= {SpecificModels._additive, SpecificModels._epistatic, 
			SpecificModels._compensatory, SpecificModels._heterogenous, SpecificModels._compound};
	public static final double logistic_odd_cutoff=0.5;
	public static final double nagative_response = -1;
	public static final double positive_response = 1;
	

	/*
	 * f(X) = \Sigma x_i
	 * Coefficients will be generated randomly using a uniform distribution 
	 */
	public static double[] additive(double[][] X) {
		int sample_size = X.length;
		int num_var = X[0].length;
		double[] response=new double[sample_size];
		double[] coef = new double[num_var];
		for(int k=0;k<num_var;k++) // randomly assign coefficients 	
			coef[k]=generator.nextGaussian();
		for(int i=0;i<sample_size;i++) {
			for(int k=0;k<num_var;k++) {
				response[i]+=coef[k]*X[i][k];
			}
		}return response;
	}
	
	/*
	 * The epistatic model, same as the AND logical operator. 
	 * 
	 * More than two vars are allowed, however suggested to just set as two (X[0].length ==2)
	 */
	public static double[] epistatic(double[][] X) {
		int sample_size = X.length;
		int num_var = X[0].length;
		double[] response=new double[sample_size];
		Arrays.fill(response, SpecificModels.positive_response);
		for(int i=0;i<sample_size;i++) {
			for(int k=0;k<num_var;k++) {
				if(X[i][k]<=0) {
					response[i]=SpecificModels.nagative_response;
					break;
				}
			}
		}return response;
	}
	
	/*
	 * The compensatory model, same as the XOR logical operator. 
	 * 
	 * Only two vars are allowed (X[0].length ==2)
	 */
	public static double[] compensatory(double[][] X) {
		int sample_size = X.length;
		int num_var = X[0].length;
		double[] response=new double[sample_size];
		if(num_var==1) {
			//Debug System.out.println("WARNING: compensatory model only supports 2 variables! But there is only 1 variable: will return the sign of this variable.");
			for(int i=0;i<sample_size;i++) {
				if((X[i][0]>0))
					response[i]=SpecificModels.positive_response;
				else response[i]=SpecificModels.nagative_response;
			}
		}else {
			if(num_var>2) {
				//Debug System.out.println("WARNING: compensatory model only supports 2 variables! The first two are used.");
			}
			for(int i=0;i<sample_size;i++) {
				if((X[i][0]>0 && X[i][1]>0)||(X[i][0]<=0 && X[i][1]<=0))
					response[i]=SpecificModels.positive_response;
				else response[i]=SpecificModels.nagative_response;
			}
		}
		return response;
	}
	
	/*
	 * The heterogenous model, same as the OR logical operator. 
	 */
	public static double[] heterogenous(double[][] X) {
		int sample_size = X.length;
		int num_var = X[0].length;
		double[] response=new double[sample_size];
		Arrays.fill(response, SpecificModels.nagative_response);
		for(int i=0;i<sample_size;i++) {
			for(int k=0;k<num_var;k++) {
				if(X[i][k]>0) {
					response[i]=SpecificModels.positive_response;
					break;
				}
			}
		}
		return response;
	}

	/*
	 * Generating the term contributed by many small (however unknown) contributors
	 * Note that double[][] grm is NOT used yet. Here use the sum of many small contributors. 
	 */
	public static double[] add_infinitesimal_term(MainFrame main_frame, double[] original, double vc) {
		double[][] vars = CausalTerm.sample_genotype_vars_wg(main_frame.num_infinitesimal, main_frame);
		double[] inf_value=new double[main_frame.num_subj_N];
		double[] coefs=new double[main_frame.num_infinitesimal];
		for(int m=0;m<main_frame.num_infinitesimal;m++) {
			coefs[m]=main_frame.generator.nextGaussian();
			for(int n=0;n<main_frame.num_subj_N;n++) {
				inf_value[n]+= coefs[m]*vars[m][n];
			}
		}
		SpecificModels.standardization_with_weight(inf_value, SpecificModels.default_weight);
		double var_ori=variance(original);
		double var_residue=var_ori*(1-vc)/vc;
		//NormalDistribution normal = new NormalDistribution(0, Math.sqrt(var_residue));
		double[] added= new double[original.length];
		for(int i=0;i<original.length;i++) 
			added[i]=original[i]+inf_value[i]*Math.sqrt(var_residue); //normal.sample();
		return added;
	} 
	
	/*
	 * Adding the residue term contributed by noise or environmental factors
	 * double h2 is the variance component of the original term; (1-vc is the variance component of the residue term.)
	 * 
	 * the variance components will ensure that var(ori)/var(total) = h2
	 * leading to the solution of var(residue)= var(ori)(1-h2)/h2.
	 * 
	 * Note that both noise and infinitesmal terms can be added using this function. 
	 */
	public static double[] add_residue_term(double[] original, double vc) {
		double var_ori=variance(original);
		double var_residue=var_ori*(1-vc)/vc;
		//NormalDistribution normal = new NormalDistribution(0, Math.sqrt(var_residue));
		double[] added= new double[original.length];
		for(int i=0;i<original.length;i++) 
			added[i]=original[i]+generator.nextGaussian()*Math.sqrt(var_residue); //normal.sample();
		return added;
	} 
	
	/*
	 * OLD VERSION using a file. Has been replaced with the Object-based.  
	 * 
	 * compound model: More complex models that are compounds of the above four models 
	 * (i.e., linear, epistatic, compensatory, and heterogenous) can also be specified. 
	 * This is achieved by recursive use of these models. For in-stance, one can first 
	 * define a “super-gene” depending on three genotypic variants using epistatic model, 
	 * and then consider this “super-gene” as one term in another function. 
	 * 
	 * Typically, the users specify the compound mode in a file formatted below:
	 * #ID	model	contributing_terms (separated by '\t')
	 * C1	additive	X1,X2,X5 
	 * C2	epistatic	C1,X3
	 * C3	compensatory	X0,X4
	 * Y	additive	C2,C3
	 * 
	 * The above forms an outcome Y=additive(epistatic(additive(X1,X2,X5),X3),compensatory(X0,X4))
	 * The final mode, Y, is the ultimate outcome. 
	 * 
	 * Note: In practice, it is recommended not to use too many terms without biological interpretations 
	 */
	public static double[] compound_file_based(double[][] X, String compound_file) { 
		int sample_size = X.length;
		int num_var = X[0].length;
		HashSet<String> supported_models=new HashSet<String>();
		for(int m=0;m<SpecificModels.supported_models.length;m++) 
			supported_models.add(SpecificModels.supported_models[m]);
		HashMap<String, double[]> all_terms=new HashMap<String, double[]>();
		for(int t=0;t<num_var;t++) {
			all_terms.put("X"+t, X[t]);
		}
		try {
			BufferedReader br=new BufferedReader(new FileReader(compound_file));
			String line=br.readLine();
			while(line.startsWith("#"))line=br.readLine(); // skip the headers
			while(line!=null) {
				String[] info=line.split("\t"); // ID, model, terms 
				String[] contributing_term_IDs=info[2].split(",");
				double[][] input_X=new double[sample_size][contributing_term_IDs.length];
				for(int i=0;i<sample_size;i++) {
					for(int k=0;k<contributing_term_IDs.length;k++) {
						input_X[i][k]=all_terms.get(contributing_term_IDs[k])[i];
					}
				}
				double[] the_new_term=response(info[1], input_X, SpecificModels.default_weight);
				all_terms.put(info[0], the_new_term);
				line=br.readLine();
			}br.close();
		}catch(Exception e) {e.printStackTrace();}
		if(all_terms.containsKey("Y")) {
			return all_terms.get("Y").clone();
		}else {
			System.out.println("Error: The last line of model file has to be Y for output.");
			return null;
		}
	} 
	
	/*
	 * compound model: More complex models that are compounds of the above four models 
	 * (i.e., linear, epistatic, compensatory, and heterogenous) can also be specified. 
	 * This is achieved by recursive use of these models. For in-stance, one can first 
	 * define a “super-gene” depending on three genotypic variants using epistatic model, 
	 * and then consider this “super-gene” as one term in another function. 
	 * 
	 * Typically, the users specify the compound mode in a file formatted below:
	 * #ID	model	contributing_terms (separated by '\t')
	 * C1	additive	X1,X2,X5 
	 * C2	epistatic	C1,X3
	 * C3	compensatory	X0,X4
	 * Y	additive	C2,C3
	 * 
	 * The above forms an outcome Y=additive(epistatic(additive(X1,X2,X5),X3),compensatory(X0,X4))
	 * The final mode, Y, is the ultimate outcome. 
	 * 
	 * Note: In practice, it is recommended not to use too many terms without biological interpretations 
	 */
	public static double[] compound(double[][] X, CompoundProtocol compound_protocol) {
		int sample_size = X.length;
		int num_var = X[0].length;
		HashSet<String> supported_models=new HashSet<String>();
		for(int m=0;m<SpecificModels.supported_models.length;m++) 
			supported_models.add(SpecificModels.supported_models[m]);
		HashMap<String, double[]> all_terms=new HashMap<String, double[]>();
		for(int t=0;t<num_var;t++) {
			double[] transpose_X_t=new double[sample_size];  // has to transpose a row into a column for X_t
			for(int i=0;i<sample_size;i++) {
				transpose_X_t[i]=X[i][t]; 
			}
			all_terms.put("X"+t, transpose_X_t);
		}
		for(int s=0;s<compound_protocol.num_steps;s++) {
			//String[] info=line.split("\t"); // ID, model, terms 
			//String[] contributing_term_IDs=compound_protocol.terms[s];
			if(compound_protocol.step_models[s].equals("compound")) {
				System.out.println("Error: There is no support for nested compound.");
				return null;
			}
			double[][] input_X=new double[sample_size][compound_protocol.step_terms[s].length];
			for(int i=0;i<sample_size;i++) {
				for(int k=0;k<compound_protocol.step_terms[s].length;k++) {
					input_X[i][k]=all_terms.get(compound_protocol.step_terms[s][k])[i];
				}
			}
			double[] the_new_term=response(compound_protocol.step_models[s], input_X, SpecificModels.default_weight);
			all_terms.put(compound_protocol.step_IDs[s], the_new_term);
		}			
		if(all_terms.containsKey("Y")) {
			return all_terms.get("Y").clone();
		}else {
			System.out.println("Error: The last line of model file has to be Y for output.");
			return null;
		}
	}
	
	/* 
	 * Randomly generate the CompoundProtocol and run it. TODO
	 */
	public static double[] compound(double[][] X) {
		 CompoundProtocol compound_protocol =new CompoundProtocol(X[0].length);
		 return compound(X, compound_protocol);
	}
	
	/*
	 * Based on models, call different functions
	 * "additive", "epistatic", "compensatory", "heterogenous", "compound"
	 */
	public static double[] response(String model, double[][] X, double weight) {
		double[] output_response;
		if(model.equals(SpecificModels._additive)) {
			output_response = SpecificModels.additive(X);
		}else if(model.equals(SpecificModels._epistatic)) {
			output_response =  SpecificModels.epistatic(X);
		}else if(model.equals(SpecificModels._compensatory)) {
			output_response =  SpecificModels.compensatory(X);
		}else if(model.equals(SpecificModels._heterogenous)) {
			output_response =  SpecificModels.heterogenous(X);
		}else if(model.equals(SpecificModels._compound)) {
			output_response =  SpecificModels.compound(X);
		}else {
			System.out.println("Error: Model "+model+" is not supported.");
			System.out.println("Only the following models are supported:\n");
			for(int m=0;m<SpecificModels.supported_models.length;m++) 
				System.out.print(SpecificModels.supported_models[m]+"; ");
			System.out.println();
			return null;
		}
		SpecificModels.standardization_with_weight(output_response, weight); // standardize to mean=0 and sd=1.
		return output_response;
	}
	
	/*
	 * given the "raw" values (containing specified biological terms only), add infinitesimal and noise terms based on
	 * predefined variance component. And finally, if specified, convert to binary according to the mode ("logistic", "liability").
	 */
	public static double[] add_inf_noise_bin(double[] response_value, MainFrame main_frame, CausalTerm causal_term) {
		double[] inf_added, noise_added;
		// add infinitesimal if it is not NaN
		if(!Double.isNaN(causal_term.infinitesimal_vc_value)) {
			double adjusted_inf_vc = causal_term.infinitesimal_vc_value/(1-causal_term.noise_vc_value);
			inf_added=SpecificModels.add_infinitesimal_term(main_frame, response_value, adjusted_inf_vc);
			noise_added = SpecificModels.add_residue_term(inf_added, causal_term.noise_vc_value);
		}
		// add noise (unconditionally)
		noise_added = SpecificModels.add_residue_term(response_value, causal_term.noise_vc_value);
		if(causal_term.ID.startsWith("Binary_Trait_")) {
			return SpecificModels.quantitative2binary(noise_added);
		}else return noise_added;		
	}
	/*
	 * Mode can be in {"logistic", "liability"}
	 */
	public static double[] quantitative2binary(double[] quantitative, String mode){
		int sample_size=quantitative.length;
		double[] binary=new double[sample_size];
		if(mode.equals("logistic")) {
			for(int i=0;i<sample_size;i++) {
				double exp=Math.pow(Math.E, quantitative[i]);
				binary[i]=(exp/(1+exp)>SpecificModels.logistic_odd_cutoff)?1:0; 
			}
		}else if(mode.equals("liability")) {
			double[] copy=quantitative.clone();
			Arrays.sort(copy);
			double cutoff=copy[sample_size/2]; // the median 
			for(int i=0;i<sample_size;i++) {
				binary[i]=(quantitative[i]>cutoff)?1:0;
			}
		}else {
			System.out.println("Error: Mode "+mode+" is invalid. Only logistic and liability are allowed.");
		}
		return binary;
	}
	
	/*
	 * Mode can be in {"logistic", "liability"}; 
	 * randomly assign a Mode based on the proportions in MainFrame.binary_mode_proportion
	 */
	public static double[] quantitative2binary(double[] quantitative){
		String mode=(SpecificModels.generator.nextDouble()<MainFrame.binary_mode_proportion[0])?
				MainFrame.binary_modes[0]:MainFrame.binary_modes[1];
		return quantitative2binary(quantitative, mode);
	}
	
	/*
	 * Calculate variance of a vector
	 */
	public static double variance(double[] data) {
		double sum = 0;
		for (int i = 0; i < data.length; i++) sum += data[i];
		double mean = sum/data.length;
		double variance = 0;
		for (int i = 0; i < data.length; i++) variance += (data[i] - mean)*(data[i] - mean);
		variance = variance/data.length;
		return variance;		
	}
	
	/*
	 * Calculate variance of a vector
	 */
	public static double correlation(double[] data1, double[] data2) {
		if(data1.length!=data2.length) {
			System.out.println("Error: correlation(data1,data2): data1.length!=data2.length");
			return Double.NaN;
		}
		int length=data1.length;
		double sum1 = 0, sum2=0;
		for (int i = 0; i < length; i++) {
			sum1 += data1[i];
			sum2 += data2[i];
		}
		double mean1 = sum1/length;
		double mean2 = sum2/length;
		double numerator = 0;
		double denominator1=0, denominator2=0;
		for (int i = 0; i < length; i++) {
			numerator += (data1[i]-mean1)*(data2[i]-mean2);
			denominator1 += (data1[i] - mean1)*(data1[i] - mean1);
			denominator2 += (data2[i] - mean2)*(data2[i] - mean2);
		}
		double correlation= numerator/(Math.sqrt(denominator1)*Math.sqrt(denominator2));
		return correlation;		
	}
	
	/*
	 * Standardization of an array (so that mean = 0, sd =1) and assign an weight
	 * 
	 * Note that if the original variance==0, return a vector with all 0s, instead of NaN. 
	 */
	public static void standardization_with_weight(double[] data, double weight) {
		double sum = 0;
		for (int i = 0; i < data.length; i++) 
			sum += data[i];
		double mean = sum/data.length;
		double variance = 0;
		for (int i = 0; i < data.length; i++) 
			variance += (data[i] - mean)*(data[i] - mean);
		variance = variance/data.length;
		if(variance==0) {
			Arrays.fill(data, 0.0);
		}else {
			for (int i = 0; i < data.length; i++) {
				data[i]=(data[i]-mean)/Math.sqrt(variance);  //Standardization to mean zero and variance 1.
				data[i]=data[i]*weight;		// multiply the weight
			}
		}		
	}
	
	/*
	 * re-scale an array's sum to be 1. 
	 */
	public static void to_prob_distribution(double[] input) {
		double sum=0;
		for(int k=0;k<input.length;k++)sum+=input[k];
		if(sum==0)System.out.println("Error: Sum of an Propotion Array = 0.");
		for(int k=0;k<input.length;k++)input[k]=input[k]/sum;
	}
	
	/*
	 * re-scale an array's sum to be 1. 
	 */
	public static double[] to_prob_distribution(int[] input) {
		double sum=0;
		for(int k=0;k<input.length;k++)sum+=input[k];
		if(sum==0)System.out.println("Error: Sum of an Propotion Array = 0.");
		double[] output=new double[input.length];
		for(int k=0;k<input.length;k++)output[k]=input[k]/sum;
		return output;
	}
}
