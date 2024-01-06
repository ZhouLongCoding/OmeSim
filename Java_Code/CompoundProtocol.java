package omesim;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;

/*
	 * compound model: More complex models that are compounds of the above four models 
	 * (i.e., additive, epistatic, compensatory, and heterogenous) can also be specified. 
	 * This is achieved by recursive use of these models. For in-stance, one can first 
	 * define a “super-gene” depending on three genotypic variants using epistatic model, 
	 * and then consider this “super-gene” as one term in another function. 
	 * 
	 * Typically, the users specify the compound mode in a file formatted below:
	 * #Step_ID	Model	Contributing_terms (separated by '\t')
	 * C1	additive	X1,X2,X5 
	 * C2	epistatic	C1,X3
	 * C3	compensatory	X0,X4
	 * Y	additive	C2,C3
	 * 
	 * The above forms a function Y=additive(epistatic(additive(X1,X2,X5),X3),compensatory(X0,X4))
	 * The final mode, Y, is the ultimate outcome. 
	 * 
	 * Note: In practice, it is recommended not to use too many terms without biological interpretations 
	 *
 */
public class CompoundProtocol {
	
	public static final String[] all_models= {"additive", "epistatic", "compensatory", "heterogenous"};
	public static final String[] nonlinear_models= {"epistatic", "compensatory", "heterogenous"};
	
	public static final int random_seed =1;
	public static final Random generator = new Random(random_seed);
	
	public int num_steps;
	public String[] step_IDs ;
	public String[] step_models ;
	public String[][] step_terms ;

	static int single_additive_cutoff=3; // using additive if num_terms is lower than simple_cutoff: num_steps = 1
	static int two_step_cutoff=5; 	// additive + 
									//one out of epiepistatic/compensatory/heterogenous: num_steps = 2
	static int three_step_cutoff=8; // additive +
									// one out of epiepistatic/compensatory/heterogenous +
									// another one out of epiepistatic/compensatory/heterogenous num_steps = 3
	//if num_terms > three_step_cutoff, then num_steps = 4
	//		additive + one out of epiepistatic/compensatory/heterogenous = C1
	//  	tow X variables using one out of epiepistatic/compensatory/heterogenous = C2
	//		Y = C1 + C2
	
	/*
	 * generating a CompoundProtocol randomly by seeding the number of X terms. 
	 * 
	 * It is not completely random. Instead the rules based on 
	 * the comments of for the three cutoffs above are followed 
	 * 
	 * No nested compound supported 
	 */
	public CompoundProtocol(int num_terms) {
		if(num_terms<2) {
			//Debug System.out.println("Warning: CompoundProtocol(int num_terms): num_terms="+num_terms+", which is < 2. Returned a single term with 'additive'.");
			this.num_steps=1;
			this.step_IDs=new String[1]; 
			this.step_IDs[0]="Y"; 
			this.step_models=new String[1]; 
			this.step_models[0]="additive"; 
			this.step_terms=new String[1][1]; 
			this.step_terms[0][0]="X0";
		}else if(num_terms <=CompoundProtocol.single_additive_cutoff) { // simply additive
			this.num_steps=1;
			this.step_IDs=new String[1]; 
			this.step_IDs[0]="Y"; 
			this.step_models=new String[1]; 
			this.step_models[0]="additive"; 
			this.step_terms=new String[1][num_terms]; 
			for(int term_index=0;term_index<num_terms;term_index++) {
				this.step_terms[0][term_index]="X"+term_index;
			}
		}else if(num_terms <=CompoundProtocol.two_step_cutoff) { // additive + epiepistatic/compensatory/heterogenous
			this.num_steps=2;
			int step2_index=CompoundProtocol.generator.nextInt(num_terms);
			this.step_IDs=new String[2]; 
			this.step_IDs[0]="C1"; 
			this.step_IDs[1]="Y"; 
			this.step_models=new String[2]; 
			this.step_models[0]="additive"; 
			this.step_models[1]=CompoundProtocol.sample_a_nonlinear(); 
			this.step_terms=new String[2][]; 
			// step 1 terms
			this.step_terms[0]=new String[num_terms-1];
			int term0_index=0;
			for(int term_index=0;term_index<num_terms;term_index++) {
				if(term_index!=step2_index)
					this.step_terms[0][term0_index++]="X"+term_index;
			}
			// step 2 terms
			this.step_terms[1]=new String[2];
			this.step_terms[1][0]="X"+step2_index;
			this.step_terms[1][1]="C1";
		}else if(num_terms <=CompoundProtocol.three_step_cutoff) {
			this.num_steps=3;
			int step2_index=CompoundProtocol.generator.nextInt(num_terms);
			int step3_index=CompoundProtocol.generator.nextInt(num_terms);
			while(step3_index==step2_index) step3_index=CompoundProtocol.generator.nextInt(num_terms);
			this.step_IDs=new String[3]; 
			this.step_IDs[0]="C1"; 
			this.step_IDs[1]="C2"; 
			this.step_IDs[2]="Y"; 
			this.step_models=new String[3]; 
			this.step_models[0]="additive"; 
			this.step_models[1]=CompoundProtocol.sample_a_nonlinear(); 
			this.step_models[2]=CompoundProtocol.sample_a_nonlinear();
			this.step_terms=new String[3][];
			// step 1 terms
			this.step_terms[0]=new String[num_terms-2];
			int term0_index=0;
			for(int term_index=0;term_index<num_terms;term_index++) {
				if(term_index!=step2_index && term_index!=step3_index)
					this.step_terms[0][term0_index++]="X"+term_index;
			}
			// step 2 terms 
			this.step_terms[1]=new String[2];
			this.step_terms[1][0]="X"+step2_index;
			this.step_terms[1][1]="C1";
			// step 3 terms
			this.step_terms[2]=new String[2];
			this.step_terms[2][0]="X"+step3_index;
			this.step_terms[2][1]="C2";
		}else {  // must be > CompoundProtocol.three_step_cutoff
			this.num_steps=4;
			int step2_index=CompoundProtocol.generator.nextInt(num_terms);
			int step3_index1=CompoundProtocol.generator.nextInt(num_terms);
			while(step3_index1==step2_index) step3_index1=CompoundProtocol.generator.nextInt(num_terms);
			int step3_index2=CompoundProtocol.generator.nextInt(num_terms);
			while(step3_index2==step2_index || step3_index2==step3_index1) 
				step3_index2=CompoundProtocol.generator.nextInt(num_terms);
			this.step_IDs=new String[4]; 
			this.step_IDs[0]="C1"; 
			this.step_IDs[1]="C2"; 
			this.step_IDs[2]="C3"; 
			this.step_IDs[3]="Y"; 
			this.step_models=new String[4]; 
			this.step_models[0]="additive"; 
			this.step_models[1]=CompoundProtocol.sample_a_nonlinear(); 
			this.step_models[2]=CompoundProtocol.sample_a_nonlinear();
			this.step_models[3]="additive"; 
			this.step_terms=new String[4][]; 
			// step 1 terms
			this.step_terms[0]=new String[num_terms-3];
			int term0_index=0;
			for(int term_index=0;term_index<num_terms;term_index++) {
				if(term_index!=step2_index && 
						term_index!=step3_index1 && term_index!=step3_index2)
					this.step_terms[0][term0_index++]="X"+term_index;
			}
			// step 2 terms
			this.step_terms[1]=new String[2];
			this.step_terms[1][0]="X"+step2_index;
			this.step_terms[1][1]="C1";
			// step 3 terms
			this.step_terms[2]=new String[2];
			this.step_terms[2][0]="X"+step3_index1;
			this.step_terms[2][1]="X"+step3_index2;
			// step 4 terms
			this.step_terms[3]=new String[2];
			this.step_terms[3][0]="C2";
			this.step_terms[3][1]="C3";
		}
		
	}
	
	public static String sample_a_nonlinear() {
		int total_num_models=CompoundProtocol.nonlinear_models.length;  // it is actually 3.
		return CompoundProtocol.nonlinear_models[CompoundProtocol.generator.nextInt(total_num_models)];
	}
	
	/*
	 * The following three methods reset the three cutoffs 
	 */
	public static void set_single_additive_cutoff(int single_additive_cutoff) {
		CompoundProtocol.single_additive_cutoff=single_additive_cutoff;
	}
	public static void set_two_step_cutoff(int two_step_cutoff) {
		CompoundProtocol.two_step_cutoff=two_step_cutoff;
	}
	public static void set_three_step_cutoff(int three_step_cutoff) {
		CompoundProtocol.three_step_cutoff=three_step_cutoff;
	}
	
	/*
	 * output the compound model to a file for future record. 
	 */
	public void out2file(String output_file) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
			bw.write("#Step_ID\tModel\tContributing_terms\n");
			for(int step_index=0;step_index<this.num_steps;step_index++) {
				bw.write(this.step_IDs[step_index]+"\t"+this.step_models[step_index]+"\t");
				for(int term_index=0;term_index<this.step_terms[step_index].length;term_index++) {
					bw.write(this.step_terms[step_index][term_index]);
					bw.write((term_index==this.step_terms[step_index].length-1)?"\n":",");
				}
			}
			bw.close();
		}catch(Exception e) {e.printStackTrace();}
	}
}
