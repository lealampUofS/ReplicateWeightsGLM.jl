module ReplicateWeightsGLM

import  GLM.stderror, 
        GLM.GlmResp, 
        GLM.LinPred, 
        GLM.AbstractGLM, 
        GLM.UnivariateDistribution, 
        GLM.Link,
        GLM.FormulaTerm,
        GLM.glm,
        GLM.coef,
        GLM.GeneralizedLinearModel,
        StatsModels.TableRegressionModel

export GeneralizedLinearModelReplicateWeights, stderror, glm_rw

mutable struct GeneralizedLinearModelReplicateWeights{G<:GlmResp,L<:LinPred} <: AbstractGLM
    rr::G
    pp::L
    fit::Bool
    maxiter::Int
    minstepfac::Float64
    atol::Float64
    rtol::Float64
    stderrors::Vector{Float64}
end

function glm_rw(formula, data; wts::Symbol, distr::UnivariateDistribution, link::Link = canonicallink(distr), numReplicatWeights::Int64=0, replicateWeightPrefix::String="", replciateWeightSuffix::String="", fayValue::Float64=0.0, kwargs...)
    main_model = glm(formula, data, distr, link; wts=data[!, wts])

    if numReplicatWeights == 0
            return main_model
    else
        mm_coef = coef(main_model)

        numCoef = length(mm_coef)

        rw_coef = Array{Vector{Float64}}(undef, 0)

        for i = 1:numReplicatWeights
            rw_symbol = Symbol(replicateWeightPrefix * string(i) * replciateWeightSuffix)
            rw_model = glm(formula, data, distr, link; wts=data[!, rw_symbol])
            push!(rw_coef, coef(rw_model))
        end

        rw_sums = zeros(numCoef)

        for coefIdx = 1:numCoef
            for i = 1:numReplicatWeights
                rw_sums[coefIdx] += (rw_coef[i][coefIdx] - mm_coef[coefIdx])^2
            end
        end
        
        se_multiplier = 1 / (numReplicatWeights * ( (1 - fayValue)^2))
        
        rw_variances = se_multiplier .* rw_sums

        stderrors = sqrt.(rw_variances)

        wm = main_model.model # wrapped main_model GLM shorthand 

        final_GLM_model = GeneralizedLinearModelReplicateWeights(wm.rr, wm.pp, wm.fit, wm.maxiter, wm.minstepfac, wm.atol, wm.rtol, stderrors)

        final_model = TableRegressionModel(final_GLM_model, main_model.mf, main_model.mm)

        return final_model
    end
end

function stderror(x::GeneralizedLinearModelReplicateWeights)
    return x.stderrors
end


end # module ReplicateWeightsGLM
