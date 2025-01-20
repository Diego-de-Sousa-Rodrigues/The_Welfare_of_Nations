The code solves the steady state in Julia and the compute the weights for France presented in the paper The Welfare of Nations: Social Preferences and the Macroeconomy.
The steady state of the model is computed using the Endogenous Grid Method (EGM), building on Alisdair McKay's code (all errors remain ours).

Julia Code for Steady-State Solution and Computation of the Weights

- Main File: `Main.jl`
  - This file is the entry point to reproduce all results presented in the paper for France.
  - Dependencies:
    - `Aiyagari_solve_endo.jl`: Solves the Aiyagari model using the Endogenous Grid Method and progressive labor taxes.
    - `parameters_France_final.jl`: Contains parameters to compute macroeconomic allocations and inequality in the steady state for France.
    - `Projection_truncation_endo.jl`: Solves the model using the truncation method.
    - `Function_pareto_weights_G_optimal.jl`: Computes the parametric and non-parametric weights.
    - `Function_Weight_Matrix.jl`: Finds individual welfare functions using political weights.
    - `Function_pareto_weights_Bin.jl`: Computes the non-parametric weights by History.
  - Additional Parameters Files:
    - `parameters_France_final_beta.jl`: Calibration for France with the U.S. discount factor.
    - `parameters_France_final_beta_taxes.jl`: Calibration for France with the U.S. discount factor and taxes.
    - `parameters_France_final_beta_inc_taxes.jl`: Calibration for France with the U.S. discount factor, taxes, and income process.
    - `parameters_France_final.jl`: Tax system of the United States applied to the French calibration to evaluate the influence of the tax system on weights.
  - Input Data:
    - `data_France.csv`: Voter turnout data for France, sourced from IPSOS data for participation rates as a function of occupations and DADS 2007 to obtain the annual income for each occupation.
  - Output:
    - Run the file `Main.jl` and ensure that all the other files are in the same folder. The output will be the figures presented in the paper for the case of France. 

Feel free to adjust this structure for your purposes or let us know if further details are needed!
