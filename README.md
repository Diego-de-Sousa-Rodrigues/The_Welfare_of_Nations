# The Welfare of Nations: Social Preferences and the Macroeconomy
**Xavier Ragot, François Le Grand, and Diego Rodrigues**

Unified Documentation for Running the Codes for France and the United States

This repository contains the codes to solve the steady state and compute the weights presented in the paper **The Welfare of Nations: Social Preferences and the Macroeconomy.** The steady state of the model is computed using the **Endogenous Grid Method (EGM)**, building on Alisdair McKay's code (all errors remain ours).

- **Structure of the Repository**
  
	- Folder:  `Codes France`
		- Contains all files necessary to reproduce the results for France.
  	- Folder: `Codes United States`
  		- Contains all files necessary to reproduce the results for the United States.
  	   
- **Common Codes**
  
The following codes are used in both Codes France and Codes United States:
	- `Main.jl`
		- The entry point to reproduce the results presented in the paper for each country. Run this file to generate the output figures. 
	- `Aiyagari_solve_endo.jl`
		- Solves the Aiyagari model using the Endogenous Grid Method and progressive labor taxes. 
	- `Projection_truncation_endo.jl'
		- Solves the model using the truncation method. 
	- `Function_pareto_weights_G_optimal.jl`
		- Computes the parametric and non-parametric weights. 
	- `Function_Weight_Matrix.jl'
		- Finds individual welfare functions using political weights. 
	- `Function_pareto_weights_Bin.jl`
		- Computes the non-parametric weights by history. 
  
- **Country-Specific Files**
  
In each folder (`Codes France` and `Codes United States`), additional files contain the calibration parameters and input data required to reproduce the results for the respective country:

	- France
		- Calibration Files:
			- `parameters_France_final.jl` 
			- `parameters_France_final_beta.jl` 
			- `parameters_France_final_beta_taxes.jl` 
			- `parameters_France_final_beta_inc_taxes.jl` 
	 	- Input Data:
			- `data_France.csv`: Voter turnout data sourced from IPSOS participation rates and DADS 2007 income data. 
  	 - United States
		- Calibration Files:
			- `parameters_US_final.jl` 
			- `parameters_US_final_beta.jl`
			- `parameters_US_final_beta_taxes.jl`
			- `parameters_US_final_beta_inc_taxes.jl`
		- Input Data:
			- `data_US.csv`: Voter turnout data sourced from Table 8 of the November 2008 Current Population Survey by the U.S. Census Bureau. 
   
- **Running the Code**
	- Place all files for the respective country in the same folder (Codes France or Codes United States). 
	- Open `Main.jl` in Julia. 
	- Execute `Main.jl` to generate the output figures presented in the paper for the respective country.

 
**Notes**
Each country folder is self-contained and includes all the necessary calibration files and data. Common codes are shared but are tailored to work with the specific parameters and data for each country.

Xavier Ragot: xavier.ragot@sciencespo.fr
François Le Grand: francois.le-grand@rennes-sb.com
Diego Rodrigues: diego.desousarodrigues@sciencespo.fr
