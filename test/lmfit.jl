import LsqFit
using Test

mutable struct Param
    name::String
    val::Float64 # will start as NaN
    min::Float64
    max::Float64
    vary::Bool
    unc::Union{Nothing, Float64}
end
function (p::Param)(;val=nothing, min=nothing, max=nothing, vary=nothing)
    isnothing(val)  || (p.val = val)
    isnothing(min)  || (p.min = min)
    isnothing(max)  || (p.max = max)
    isnothing(vary) || (p.vary = vary)
    return nothing
end
uninit_param(name::String) = Param(name, NaN, -Inf, Inf, true, nothing)

function Base.getindex(q::Vector{Param}, key::String) #emulate an ordered dict
    for v in q
        v.name == key && return v
    end
    KeyError("$key not among possible keys $([v.name for v in q])")
end
# not type stable, but its just for the repl so shouldnt matter too much
# this is type piracy, for a quick and dirty demo
# should make a custom container and implement the array interface
Base.propertynames(q::Vector{Param}) = tuple([Symbol(p.name) for p in q]...)
Base.getproperty(q::Vector{Param}, key::Symbol) = q[string(key)]
Base.values(q::Vector{Param}) = [p.val for p in q]
Base.keys(q::Vector{Param}) = [p.name for p in q]

function Base.show(io::IO, p::Param)
    unc = isnothing(p.unc) ? "NA" : p.unc
    print(io, "Param(\"$(p.name)\", val=$(p.val), min=$(p.min), max=$(p.max), vary=$(p.vary), unc=$(unc))")
end

struct Result
    model
    x::Vector{Float64}
    init_params::Vector{Param}
    params::Vector{Param}
    fit::LsqFit.LsqFitResult
    weights::Union{Nothing, Vector{Float64}}
end
(r::Result)(;x=r.x, params=r.params) =  r.model(x=x, params=params)

# add parameter guessing
struct Model
    name::String
    f
    n_params::Int
end
(m::Model)(;x, params) = m.f(x, values(params)...)
function make_model(f;name=nothing)
    params = params_from_function(f) # bit redundant with model_and_params but this shouldn't be in a hot loop
    if isnothing(name)
        return Model(string(f), f, length(params))
    else
        return Model(name, f, length(params))
    end
end
function model_and_params(f;name=nothing)
    params = params_from_function(f)
    return make_model(f,name=name), params
end
function get_reduced_fit_function(model::Model, params)
        (x, p_vary) -> begin
            p = expand_params(p_vary, params)
            model.f(x, p...)
        end    
end


# sketch of a composite model
struct CompositeModel
    binary_op
    m1::Model
    m2::Model
end
(m::CompositeModel)(;x, params) = m.binary_op(m.m1(x=x, params=params), m.m2(x=x, params=params))
function get_reduced_fit_function(m::CompositeModel, params)
    (x, p_vary) -> begin
        p = expand_params(p_vary, params)
        a = m.m1.f(x, p[1:m.m1.n_params]...)
        b = m.m2.f(x, p[m.m1.n_params+1:end])
        m.binary_op(a,b)
    end
end

function params_from_function(f)
    ms = methods(f)
    @assert length(ms) == 1 "not implemented for functions with more than one method, call params_from_method"
    params_from_method(first(ms))
end

function params_from_method(m)
    argnames = Base.method_argnames(m)
    # the first argname is self, which we don't want
    @assert argnames[2] == :x "first argument of m must be x, m=$m"
    @assert allunique(argnames) "no duplicate param names allowed"
    [uninit_param(string(name)) for name in argnames[3:end]]
end

function expand_params(p_vary, params)
    # p_vary is a vector of varying parameters values with length equal to the number of varying params
    # use params to add non-varying parameter values to create 
    # p_all, with length equal to params
    # p_all is a vector of parameter values passed to the underlying fitter
    p_all = zeros(length(params))
    i_vary = 0 # index into p_vary
    for (j, param) in enumerate(params)
        if param.vary
            i_vary+= 1
            p_all[j] = p_vary[i_vary]
        else
            p_all[j] = param.val
        end
    end
    n_vary = i_vary
    @assert n_vary == length(p_vary) "length(p_vary) = $(length(p_vary)) should be $(n_vary)"
    return p_all
end

function result_from_fit(model, x, init_params, fit, weights)
    params = deepcopy(init_params)
    sigma = LsqFit.stderror(fit)
    i_vary = 0
    for p in params
        if p.vary
            i_vary += 1
            p.val = fit.param[i_vary]
            p.unc = sigma[i_vary]
        end
    end
    Result(model, x, init_params, params, fit, weights)
end
    
function fit(model, params; x, y, weights=nothing)
    validate_params(params)
    reduced_fit_func = get_reduced_fit_function(model, params)
    p0 = [p.val for p in params if p.vary]
    lower = [p.min for p in params if p.vary]
    upper = [p.max for p in params if p.vary]
    if weights === nothing
        fit = LsqFit.curve_fit(reduced_fit_func, x, y, p0, lower=lower, upper=upper) 
    else
        fit = LsqFit.curve_fit(reduced_fit_func, x, y, weights, p0, lower=lower, upper=upper) 
    end
    result_from_fit(model, x, params, fit, weights)
end

function validate_params(params)
    hasnan(p::Param) = isnan(p.val) || isnan(p.min) || isnan(p.max)
    invalid = [p for p in params if hasnan(p)]
    if length(invalid) > 0
        error("these parameters have NaNs: $([p.name for p in invalid])")
    end
    return nothing
end



if true
# test with
x = 1:100;
ydata = x.^2.5;
f(x, a, b) = x.^a.+b
model, params = model_and_params(f)
params.b(val = 0, vary = false)
params.a(val = 3)
@test [2.5, 0] == expand_params([2.5], params)
@test model(x=x,params=params) == f(x,params.a.val,params.b.val)
reduced_fit_func = get_reduced_fit_function(model, params)
@test reduced_fit_func(x, [2.5]) == ydata # takes a vector of floats for which param varies

result = fit(model, params; x=x, y=ydata)

end