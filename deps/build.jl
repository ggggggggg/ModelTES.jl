println("running ModelTES build.jl: "*@__FILE__)

using Pkg
println("Installing ARMA.jl")
Pkg.add("https://github.com/joefowler/ARMA.jl")
Pkg.activate(".")
Pkg.build()
Pkg.test(; coverage=true)
