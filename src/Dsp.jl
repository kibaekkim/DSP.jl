# __precompile__()

module Dsp

include("DspCInterface.jl")
include("DspSolverInterface.jl")
include("block.jl")

using Dsp.DspSolverInterface
# export DspSolver
# import JuMP

end # module