# __precompile__()

module Dsp

include("DspCInterface.jl")
include("block.jl")

using Dsp.DspCInterface

import JuMP
export readSmps, getblocksolution, optimize, getprimobjval, getdualobjval, getsolutiontime, getblockids

# DspModel placeholder
model = DspModel()

###############################################################################
# JuMP.solvehook
###############################################################################

# This function is hooked by JuMP (see block.jl)
function dsp_solve(m::JuMP.Model; suppress_warnings = false, comm = nothing, options...)

    # free DspModel
    DspCInterface.freeModel(Dsp.model)

    # parse options
    for (optname, optval) in options
        if optname == :param
            DspCInterface.readParamFile(Dsp.model, optval)
        elseif optname == :solve_type
            if optval in [:Dual, :Benders, :Extensive]
                Dsp.model.solve_type = optval
            else
                warn("solve_type $optval is not available.")
            end
        else
            warn("Options $optname is not available.")
        end
    end

    # load JuMP model to Dsp
    DspCInterface.loadProblem(Dsp.model, m)

    # solve
    DspCInterface.solve(Dsp.model, comm)

    # solution status
    statcode = DspCInterface.getStatus(Dsp.model)
    stat = parseStatusCode(statcode)

    # Extract solution from the solver
    Dsp.model.numRows = DspCInterface.getTotalNumRows(Dsp.model)
    Dsp.model.numCols = DspCInterface.getTotalNumCols(Dsp.model)
    m.objVal = NaN
    m.colVal = fill(NaN, Dsp.model.numCols)

    if stat != :Optimal
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end

    if !(stat == :Infeasible || stat == :Unbounded)
        try
            getDspSolution(suppress_warnings)

            # Don't corrupt the answers if one of the above two calls fails
            m.objVal = Dsp.model.primVal
            # parse solution to each block
            n_start = 1
            n_end = m.numCols
            m.colVal = copy(Dsp.model.colVal[n_start:n_end])
            n_start += m.numCols
            if haskey(m.ext, :DspBlocks) == true
                blocks = m.ext[:DspBlocks].children
                for s in keys(blocks)
                    n_end += blocks[s].numCols
                    blocks[s].colVal = Dsp.model.colVal[n_start:n_end]
                    n_start += blocks[s].numCols
                end
            end
        end
    end

    # maximization?
    if m.objSense == :Max
        m.objVal *= -1
        dsp.primVal *= -1
        dsp.dualVal *= -1
    end

    # Return the solve status
    stat

end

###############################################################################
# Input/output files
###############################################################################

function optimize(;suppress_warnings = false, comm = nothing, options...)

    # parse options
    for (optname, optval) in options
        if optname == :param
            DspCInterface.readParamFile(Dsp.model, optval)
        elseif optname == :solve_type
            if optval in [:Dual, :Benders, :Extensive]
                Dsp.model.solve_type = optval
            else
                warn("solve_type $optval is not available.")
            end
        else
            warn("Options $optname is not available.")
        end
    end

    # solve
    DspCInterface.solve(Dsp.model, comm)

    # solution status
    statcode = DspCInterface.getStatus(Dsp.model)
    stat = parseStatusCode(statcode)

    # Extract solution from the solver
    Dsp.model.numRows = DspCInterface.getTotalNumRows(Dsp.model)
    Dsp.model.numCols = DspCInterface.getTotalNumCols(Dsp.model)

    if stat != :Optimal
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end

    if !(stat == :Infeasible || stat == :Unbounded)
        getDspSolution(suppress_warnings)
    end

    # Return the solve status
    stat
end

# Read model from SMPS files
function readSmps(filename::AbstractString)

    # free DspModel
    DspCInterface.freeModel(Dsp.model)

    # read Smps file
    DspCInterface.readSmps(Dsp.model, filename)
end

###############################################################################
# Helpers
###############################################################################

# get solutions for blocks
function getblocksolution(m::JuMP.Model)
    sol = []
    if haskey(m.ext, :DspBlocks) == true
        blocks = m.ext[:DspBlocks].children
        for s in keys(blocks)
            push!(sol, blocks[s].colVal)
        end
    else
        warn("No block was created in the model.")
    end
    sol
end

# get dual objective value
getprimobjval() = Dsp.model.primVal
# get dual objective value
getdualobjval() = Dsp.model.dualVal
# get dual value
JuMP.getdual() = Dsp.model.rowVal
# get solution time
getsolutiontime() = DspCInterface.getWallTime(Dsp.model)
# get block ids
getblockids(nblocks::Integer) = DspCInterface.getProcIdxSet(nblocks)

function parseStatusCode(statcode::Integer)
    stat = :NotSolved
    if statcode == 3000
        stat = :Optimal
    elseif statcode == 3001
        stat = :Infeasible
    elseif statcode == 3002
        stat = :Unbounded
    elseif statcode in [3004,3007,3010]
        stat = :IterOrTimeLimit
    elseif statcode == 3005
        stat = :GapTolerance
    elseif statcode == 3006
        stat = :NodeLimit
    elseif statcode in [3008,3009,3014,3015,3016]
        stat = :UserLimit
    elseif statcode in [3011,3012,3013,3999]
        stat = :Error
    else
        stat = :Unknown
        warn("Unknown status: $statcode")
    end

    stat
end

function getDspSolution(suppress_warnings)
    Dsp.model.primVal = DspCInterface.getPrimalBound(Dsp.model)
    Dsp.model.dualVal = DspCInterface.getDualBound(Dsp.model)

    if Dsp.model.solve_type == :Dual
        Dsp.model.rowVal = DspCInterface.getDualSolution(Dsp.model)
    else
        Dsp.model.colVal = DspCInterface.getSolution(Dsp.model)
    end
end

end # module