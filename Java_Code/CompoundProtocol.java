package omesim;

import java.util.ArrayList;

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
	 * The above forms a function Y=additive(epistatic(additive(X1,X2,X5),X3),compensatory(X0,X4))
	 * The final mode, Y, is the ultimate outcome. 
	 * 
	 * Note: In practice, it is recommended not to use too many terms without biological interpretations 
	 *
 */
public class CompoundProtocol {
	
	public int num_steps;
	public String[] step_IDs ;
	public String[] step_models ;
	public String[][] terms ;
	
	/*
	 * generating a CompoundProtocol randomly by seeding the number of X terms.
	 * 
	 * No nested compound supported 
	 */
	public CompoundProtocol(int num_terms) {
		
	}
	
	public void out2file(String output_file) {
		
	}
}
