# __precompile__()

module Dsp

include("DspCInterface.jl")
include("block.jl")

using Dsp.DspCInterface

import JuMP
export getdsp,
  # solve,
  readSmps

# DspModel placeholder
dsp = nothing

# return DspModel object
function getdsp()
  Dsp.dsp
end

###############################################################################
# JuMP.solvehook
###############################################################################

# This function is hooked by JuMP (see block.jl)
function dsp_solve(m::JuMP.Model; suppress_warmings = false, options...)

    dsp = Dsp.dsp
    if dsp == nothing
        # create DspModel
        dsp = DspModel()
    else
        # free DspModel
        DspCInterface.freeModel(dsp)
    end

    # parse options
    for (optname, optval) in options
        if optname == :param
            DspCInterface.readParamFile(dsp, optval)
        elseif optname == :solve_type
            if optval in [:Dual, :Benders, :Extensive]
                dsp.solve_type = optval
            else
                warn("solve_type $optval is not available.")
            end
        else
            warn("Options $optname is not available.")
        end
    end

    # load JuMP model to Dsp
    DspCInterface.loadProblem(dsp, m)

    # solve
    DspCInterface.solve(dsp)

    # solution status
    statcode = DspCInterface.getSolutionStatus(dsp)
    stat = :Error
    if statcode == 3000
        stat = :Optimal
    elseif statcode == 3001
        stat = :Infeasible
    elseif statcode == 3002
        stat = :Unbounded
    elseif statcode in [3004,3005,3006,3007,3008,3009,3010,3012,3014,3015]
        stat = :UserLimit
    elseif statcode in [3011,3013,3999]
        stat = :Error
    else
        error("Internal library error")
    end

    # Extract solution from the solver
    dsp.numRows = DspCInterface.getTotalNumRows(dsp)
    dsp.numCols = DspCInterface.getTotalNumCols(dsp)
    m.objVal = NaN
    m.colVal = fill(NaN, dsp.numCols)

    if stat != :Optimal
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end

    if !(stat == :Infeasible || stat == :Unbounded)
        try
            dsp.objVal = DspCInterface.getObjValue(dsp)
            dsp.colVal = DspCInterface.getSolution(dsp)

            if dsp.solve_type == :Dual
                # dualVal
                warn("Not implemented: getting dual solution")
            end

            # Don't corrupt the answers if one of the above two calls fails
            m.objVal = dsp.objVal
            # parse solution to each block
            n_start = 1
            n_end = m.numCols
            m.colVal = copy(dsp.colVal[n_start:n_end])
            n_start += m.numCols
            if haskey(m.ext, :DspBlocks) == true
                blocks = m.ext[:DspBlocks].children
                for s in keys(blocks)
                    n_end += blocks[s].numCols
                    blocks[s].colVal = dsp.colVal[n_start:n_end]
                    n_start += blocks[s].numCols
                end
            end
        end
    end

    # maximization?
    if m.objSense == :Max
        m.objVal *= -1
        dsp.objVal *= -1
    end

    # Return the solve status
    stat

end

###############################################################################
# Input/output files
###############################################################################

function JuMP.solve(;suppress_warmings = false, options...)

    # parse options
    for (optname, optval) in options
        if optname == :param
            DspCInterface.readParamFile(dsp, optval)
        elseif optname == :solve_type
            if optval in [:Dual, :Benders, :Extensive]
                dsp.solve_type = optval
            else
                warn("solve_type $optval is not available.")
            end
        else
            warn("Options $optname is not available.")
        end
    end

    # solve
    DspCInterface.solve(dsp)

    # solution status
    statcode = DspCInterface.getSolutionStatus(dsp)
    stat = :Error
    if statcode == 3000
        stat = :Optimal
    elseif statcode == 3001
        stat = :Infeasible
    elseif statcode == 3002
        stat = :Unbounded
    elseif statcode in [3004,3005,3006,3007,3008,3009,3010,3012,3014,3015]
        stat = :UserLimit
    elseif statcode in [3011,3013,3999]
        stat = :Error
    else
        error("Internal library error")
    end

    # Extract solution from the solver
    dsp.numRows = DspCInterface.getTotalNumRows(dsp)
    dsp.numCols = DspCInterface.getTotalNumCols(dsp)

    if stat != :Optimal
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end

    if !(stat == :Infeasible || stat == :Unbounded)
        try
            dsp.objVal = DspCInterface.getObjValue(dsp)
            dsp.colVal = DspCInterface.getSolution(dsp)

            if dsp.solve_type == :Dual
                # dualVal
                warn("Not implemented: getting dual solution")
            end
        end
    end

    # Return the solve status
    stat
end

# Read model from SMPS files
function readSmps(filename::AbstractString)

    dsp = Dsp.dsp
    if dsp == nothing
        # create DspModel
        dsp = DspModel()
    else
        # free DspModel
        DspCInterface.freeModel(dsp)
    end

    # read Smps file
    DspCInterface.readSmps(dsp, filename)
end

end # module