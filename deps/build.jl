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
  cd(Pkg.dir("ARMA"))
  ;git fetch --tags
  ;git checkout v0.1.1
end
