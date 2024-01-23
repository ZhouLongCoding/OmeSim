User’s Manual of OmeSim 1.0

Zhou Long & Qingrun Zhang
Jan 21, 2024


1.	Introduction 

OmeSim is a comprehensive simulator that simulates molecular omics data (i.e., whole-transcriptome data in its first release) and multiple (possibly correlated) phenotypic traits simultaneously. The immediate application of OmeSim is to support the assessment of novel statistical models and computational tools that aim to integrate multi-omics in the discovery of genetic basis of complex traits.  

Unlike most other simulators, OmeSim emphasizes the following features that are becoming the focus of many statistical and computational tools:

•	Simulates the whole-transcriptome and the whole-phenome together based on user-specified genotype file, gene file, pathway file, and many parameters. It supports nonlinear genetic models including epistasis, compensatory, heterogenous and the compound combination of them. This is in addition to the standard additive model and the infinitesimal model. 
•	Supports various causality models (causality, pleiotropy, and reverse causality) at the level of individual terms (i.e., genetic variants, gene expressions, and traits), forming a whole-transcriptome and whole-phenome causality graph (e.g., Figure ? in Section Output Files)  
•	Outputs expression-expression, expression-trait and trait-trait correlations as well as genetics-expression and genetics-trait associations BEFORE adding noise (= random residual) and in infinitesimal (= contribution from the genetic background) terms. 

The genuine correlation and associations, together with the causality graph, serve as “gold-standard” to assess statistical properties of tools discovering genetic basis in the presence of complicated linear and nonlinear relationship. 

In addition to the function Simulate that simulates data, another function Causality provide functions to check causality relationship between any three terms under investigation based on the simulated causality graph.  

2.	A quick start

After downloading the tar ball OmeSim.tar.gz, please just decompress the file by 

> tar -xvzf OmeSim.tar.gz 

One will see the executable OmeSim.jar, a parameter file (parameter.txt), sample input files, and sample output files. Please modify the paths of the sample input files and the output file folder corresponding to the project folder in your local computer, and type the following commend:

> java -Xmx4g -jar OmeSim.jar Simulate -input parameter.txt 
 
Then one will see the outcome files in the specified folder. By comparing them to the example outcome files, one can verify whether the program is running smoothly.  
![image](https://github.com/ZhouLongCoding/OmeSim/assets/96537327/de1d353f-3ac0-4f88-bfd2-4d21cc4b659d)
