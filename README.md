# ReplicateWeightsGLM.jl
 An extension package for GLM.jl to support analysis using replicate weights

## Replicate weight algorithms supported

:x: Bootstrap
:x: Jackknife
:white_check_mark: BRR
:white_check_mark: Fay's BRR 

## Usage
To use replicate weights for variance estimation the `glm_rw` function is used
The data is passed to the function in a single `DataFrame`. 
The replicate weight columns are specified using a suffix and prefix string. 
The replicate weight columns are assumed to be in asscending numeric order.  
The forumla is specified using a [StatsModel.jl](https://juliastats.org/StatsModels.jl/stable) [`Formula Object`](https://juliastats.org/StatsModels.jl/stable/formula/).
Most functions follow the `glm` function call specified in [GLM.jl](https://juliastats.org/GLM.jl/stable/#Fitting-GLM-models).
Below is an example using Fay's BRR.

```julia
#Note: Example code. This does not run. 
#specific function parameters need to be specified based on the data
model_formula = @formula(y ~ x)
output_model = glm_rw(model_formula, df, :wts, Binomial(), LogitLink(), 500, "rw", "", 0.5)
println(coeftable(output_model))
```