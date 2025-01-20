The code solves the steady state in Julia and the compute the weights for the United States presented in the paper The Welfare of Nations: Social Preferences and the Macroeconomy.
The steady state of the model is computed using the Endogenous Grid Method (EGM), building on Alisdair McKay's code (all errors remain ours).

Julia Code for Steady-State Solution and Computation of the Weights

- Main File: `Main.jl`
  - This file is the entry point to reproduce all results presented in the paper for the United States.
  - Dependencies:
    - `Aiyagari_solve_endo.jl`: Solves the Aiyagari model using the Endogenous Grid Method and progressive labor taxes.
    - `parameters_US_final.jl`: Contains parameters to compute macroeconomic allocations and inequality in the steady state for the United States.
    - `Projection_truncation_endo.jl`: Solves the model using the truncation method.
    - `Function_pareto_weights_G_optimal.jl`: Computes the parametric and non-parametric weights.
    - `Function_Weight_Matrix.jl`: Finds individual welfare functions using political weights.
    - `Function_pareto_weights_Bin.jl`: Computes the non-parametric weights by History.
  - Additional Parameters Files:
    - `parameters_US_final_beta.jl`: Calibration for the U.S. with the French discount factor.
    - `parameters_US_final_beta_taxes.jl`: Calibration for the U.S. with the French discount factor and taxes.
    - `parameters_US_final_beta_inc_taxes.jl`: Calibration for the U.S. with the French discount factor, taxes, and income process.
    - `parameters_France_final.jl`: Tax system of France applied to the U.S. calibration to evaluate the influence of the tax system on weights.
  - Input Data:
    - `data_US.csv`: Voter turnout data for the United States, sourced from Table 8 of the November 2008 Current Population Survey by the U.S. Census Bureau.
  - Output:
    - Run the file `Main.jl` and ensure that all the other files are in the same folder. The output will be the figures presented in the paper for the case of the United States. 

Feel free to adjust this structure for your purposes or let us know if further details are needed!
