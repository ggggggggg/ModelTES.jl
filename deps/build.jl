println("running ModelTES build.jl: "*@__FILE__)

"installed(s::String)
If the package whose name is `s` is installed return `true`, otherwise return `false`."
function installed(s::String)
  try
    a=Pkg.installed(s)
    if a == nothing
      return false
    else
      return true
    end
  catch
    return false
  end
end

if !installed("ARMA")
  println("Installing ARMA.jl")
  Pkg.clone("https://github.com/joefowler/ARMA.jl")
  # The above is the master branch of ARMA, which is no good for Julia 0.6.
  # We'll need to checkout a valid ARMA version, and then re-build this package
  cd(Pkg.dir("ARMA"))
  run(`git fetch --tags`)
  run(`git checkout v0.1.1`)
  println("Re-building ModelTES:")
  Pkg.build("ModelTES")
end
