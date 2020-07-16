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
end
(r::Result)(;x=r.x, params=r.params) =  r.model(x=x, params=params)

struct Model
    name::String
    f
end
(m::Model)(;x, params) = m.f(x, values(params)...)

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

function get_reduced_fit_func_function(f, params)
    (x, p_vary) -> begin
        p = expand_params(p_vary, params)
        f(x, p...)
    end    
end

function result_from_fit(f, x, init_params, fit)
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
    Result(f, x, init_params, params, fit)
end
    
function fit(f, params; x, y)
    validate_params(params)
    reduced_fit_func = get_reduced_fit_func_function(f, params)
    p0 = [p.val for p in params if p.vary]
    lower = [p.min for p in params if p.vary]
    upper = [p.max for p in params if p.vary]
    fit = LsqFit.curve_fit(reduced_fit_func, x, y, p0, lower=lower, upper=upper)   
    result_from_fit(f, x, params, fit)
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
params = params_from_function(f)
params.b(val = 0, vary = false)
params.a(val = 3)
@test [2.5, 0] == expand_params([2.5], params)
# result = fit(f, params, data)
reduced_fit_func = get_reduced_fit_func_function(f, params)
@test ydata == reduced_fit_func(x, [2.5])

result = fit(f, params; x=x, y=ydata)

end