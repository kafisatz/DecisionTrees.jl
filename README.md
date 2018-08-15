# DecisionTrees.jl
Decision Tree Algorithms for Regression

This package is under development

Linux: [![Build Status](https://travis-ci.org/kafisatz/DecisionTrees.jl.svg?branch=master)](https://travis-ci.org/kafisatz/DecisionTrees.jl)

Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/kafisatz/DecisionTrees.jl?branch=master&svg=true)](https://ci.appveyor.com/project/kafisatz/DecisionTrees-jl/branch/master)

[![Coverage Status](https://coveralls.io/repos/kafisatz/DecisionTrees.jl/badge.svg?branch=master)](https://coveralls.io/r/kafisatz/DecisionTrees.jl?branch=master)
[![codecov.io](http://codecov.io/github/kafisatz/DecisionTrees.jl/coverage.svg?branch=master)](http://codecov.io/github/kafisatz/DecisionTrees.jl?branch=master)


## Installation

This package only works with Julia 0.7 or higher. 
To install DecisionTrees.jl, run 
`using Pkg; pkg"add https://github.com/kafisatz/DecisionTrees.jl"` 
at the REPL. Alternatively you can enter package mode by typing `]`. In package mode you can simply run `add https://github.com/kafisatz/DecisionTrees.jl`.

We note that this package uses PyCall to create Excel files of the model results. Therefore, there are a number of dependencies (PyCall, Conda, Python). These should be installed automatically when you use the package for the first time. It might be helpful to run `pkg"test DecisionTrees"` after the installation to see whether it was installed properly. You may want to run `Pkg.build("Conda");Pkg.build("PyCall");` if either of these packages is not built properly.

## Documentation

This is in the works.

This package provides Regression Decision Trees. Classification is not possible (except for a two class problem which can be viewed as binary response variable 0/1).
Currently the following model types are provided:
* Single Decision Tree
* Boosted Decision Trees

An approach for bagging (and random forest) is in development.

We note that the entire code is for the construction of the decision trees is written in Julia.

### Motivation
This package was developed for uses cases in insurance where other decision tree packages may not always be appropriate. Specifically the response variable for the trees is always considered to be a **ratio**. That is to say the "response variable" in these models always consists of **two** variables: a numerator and a denominator. As an example consider the task of **claim frequency modelling**: If we want to predict the number of insurance claims for any given policyholder during the next year (or an arbitrary period), we need to consider two elements: the number of claims that may occur and the time that the policy holder is insured (the exposure). Both of these elements are typically present in insurance data. The modelling of **loss ratios** (i.e. the sum of losses divided by the sum of premiums) is a similar problem. Both loss ratio and claim frequency modelling are difficult to be modelled with the mean squared error as a splitting criterion (which is the default in many decision tree implementations). One might be tempted to simply evaluate the claim frequency for each row in the data set and then use it as the target variable (in any decision tree or other model), however due to the nature of the problem, this may not necessarily yield good results. 

Because the package is modelling ratios (instead of a single real valued variable)

There are likely use cases in other areas, which have a similar problem setting.

## Ideas for future improvements

* Review splitting criteria (especially poisson,gamma and mean squared error). Possibly compare them to other implementations.
* Review the formulae & calculations of the Gini and RSS metrics
* Each node should not only store 'subset' but also the actual value (being a set of Strings or a or Float64 value) value which defines the split
* Implement Bagging
* Implement getindex for DTMTable (such that dtmtable[1:30] returns another valid dtmtable struct)
* fix/add/improve certain error statistics; mainly the RSS and Gini
* Clean up the naming convention of the field names of ModelSettings (bool_xy, BoolXY, ...)
* Revise and DataStreams are in the Project.toml file but not specifically needed for the package. However we are having issues with the package manager in 0.7 alpha and beta (for now). Having the two packages in the Project file solves these issues.
* Pruning (manual pruning: add a function prune_leaf(model,leafnrs::Vector{Int}); automated pruning based on some performance metric (on trn or val data))
* consider threading and simd to improve performance

Also: 
* To fix this: the code should not rely on dtmtable.mappings and dtmtable.candmatwomax any more (if possible). The pools have the equivalent information
* Add more unit tests to trigger errors and warnings (see also the Julia docs on testing)
* Add meaningful unit tests, e.g. check the lift on trn and val for a given single tree and boosting (for a given seed)
* To add purity improvement for each split! Stop splitting if it is > some_specified_threshold
* Possibly add variable importance or partial dependence plot
* consider SSEPointwise

* Possibly add trnidx and validx to the resulting ensemble. This is relevant in case of a CV sampling which is performed. Otherwise it is not possible to reconstruct the Excel statistics after the model has run.
* Add 'time_finished' to model result (and possibly the time needed for the modelling)
* Need to ensure that the code runs smoothly even when no split is found (e.g. minw too big)
* Possibly introduce an option not to calculate certain intermediate metrics (scoring, sorting, ...) -> faster runtime
* To check if sortperm! performance has improved in newer Julia versions
* Do not use df_name_vector any more, but the symbol in the subset or simply names(features)
* Add yes/no to the top of each printed DOT graph (in the PDF output) indicating that 'yes goes left'
* Try to avoid branching in critical functions; use ifelse or similar constructs (this might improve the performance)
* Clean up code (too many type annotations, ... )
* Consider a nicer representation (such as https://github.com/sisl/D3Trees.jl)
* 'Remove' variables which are not used by a model from the generated SAS/VBA/C#code (e.g. 'dim XY' in VBA)