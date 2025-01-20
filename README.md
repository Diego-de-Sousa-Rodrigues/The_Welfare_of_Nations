# The_Welfare_of_Nations

Unified Documentation for Running the Codes for France and the United States
This repository contains the codes to solve the steady state and compute the weights presented in the paper The Welfare of Nations: Social Preferences and the Macroeconomy. The steady state of the model is computed using the Endogenous Grid Method (EGM), building on Alisdair McKay's code (all errors remain ours).
Structure of the Repository
- Folder:  `Codes France’
	- Contains all files necessary to reproduce the results for France. 
- Folder: `Codes United States’
		- Contains all files necessary to reproduce the results for the United States. 
Common Codes
The following codes are used in both Codes France and Codes United States:
- `Main.jl’
	- The entry point to reproduce the results presented in the paper for each country. Run this file to generate the output figures. 
- `Aiyagari_solve_endo.jl’
	- Solves the Aiyagari model using the Endogenous Grid Method and progressive labor taxes. 
- `Projection_truncation_endo.jl’
	- Solves the model using the truncation method. 
- `Function_pareto_weights_G_optimal.jl’
- Computes the parametric and non-parametric weights. 
- `Function_Weight_Matrix.jl’
	- Finds individual welfare functions using political weights. 
- `Function_pareto_weights_Bin.jl’
	- Computes the non-parametric weights by history.
