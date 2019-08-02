# Two-way_effect_modifications

Example R codes to assess two-way effect modifications of air pollution and air temperature on daily mortality: (1) Air pollution effects stratified by air temperature; (2) Air temperature effects stratified by air pollution concentrations. Poisson additive models with over-dispersion, together with a penalized distributed lag nonlinear temperature term developed by Dr. Antonio Gasparrini (https://github.com/gasparrini), is applied as the basic confounder model.

To examine effect modification by air temperature, air temperature is categorized into three levels. An interaction term between air pollutant and categorized air temperature at the same lag structure is added to the basic confounder model. Please use the "Modbytemp_3cat.R" function to get the air pollution effects stratified by three air temperature levels.

To examine effect modification by air pollutants, air pollutant is categorized into two levels. An interaction term between the penalized distributed lag nonlinear temperature term and an air pollutant strata indicator was added to the basic confounder model. Based on Dr. Gasparrini's DLNM package, a modified function in "crossreduce_int_2APcats.R" is applied to get air temperature effects stratified by air pollution concentrations.

For more details, please refer to the 2018 publication:

Chen K, Wolf K, Breitner S, Gasparrini A, Stafoggia M, Samoli E, et al. Two-way effect modifications of air pollution and air temperature on total natural and cardiovascular mortality in eight European urban areas. Environment International. 2018;116:186-96. doi: https://doi.org/10.1016/j.envint.2018.04.021.
