# DecisionTrees.jl
Julia Decision Tree Algorithms for Regression

This package is under development



Linux: [![Build Status](https://travis-ci.org/kafisatz/DecisionTrees.jl.svg?branch=master)](https://travis-ci.org/kafisatz/DecisionTrees.jl)

Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/kafisatz/DecisionTrees.jl?branch=master&svg=true)](https://ci.appveyor.com/project/kafisatz/DecisionTrees-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/kafisatz/DecisionTrees.jl/badge.svg?branch=master)](https://coveralls.io/r/kafisatz/DecisionTrees.jl?branch=master)
[![codecov.io](http://codecov.io/github/kafisatz/DecisionTrees.jl/coverage.svg?branch=master)](http://codecov.io/github/kafisatz/DecisionTrees.jl?branch=master)




## Ideas for future improvements

* Add MSE (normal deviance) as a splitting criterion
* Implement Bagging
* fix/add/improve certain error statistics; mainly the RSS and Gini

Also: 
* To fix this: the code should not rely on dtmtable.mappings and dtmtable.candmatwomax anymore (if possible). The pools have the equivalent information
* Add more unit tests to trigger errors and warnings (see also the Julia docs on testing)
* Add meaningful unit tests, e.g. check the lift on trn and val for a given single tree and boosting (for a given seed)
* To add purity improvement for each split! Stop splitting if it is > some_specified_threshold
* Possibly add variable importance or partial dependence plot

* Possibly add trnidx and validx to the resulting ensemble. This is relevant in case of a CV sampling which is performed. Otherwise it is not possible to reconstruct the Excel statistics after the model has run.
* Add 'time_finished' to model result (and possibly the time needed for the modelling)
* Need to ensure that the code runs smoothly even when no split is found (e.g. minw too big)
* Possibly introduce an option not to calculate certain intermediate metrics (scoring, sorting, ...) -> faster runtime
* To check if sortperm! performance has improved in newer Julia versions
* Do not use df_name_vector anymore, but the symbol in the subset or simply names(features)
* Add yes/no to the top of each printed DOT graph (in the PDF output) indicating that 'yes goes left'
* Try to avoid branching in critical functions; use ifelse or similar constructs (this might improve the performance)
* Clean up code (too many type annotations, ... )
* Consider a nicer representation (such as https://github.com/sisl/D3Trees.jl)
* 'Remove' variables which are not used by a model from the generated SAS/VBA/C#code (e.g. 'dim XY' in VBA)
