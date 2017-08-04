# __precompile__()

module Dsp

include("DspCInterface.jl")
include("block.jl")

using Dsp.DspCInterface

import JuMP
export
    readSmps, 
    getblocksolution, 
    optimize, 
    getprimobjval, 
    getdualobjval, 
    getprimvalue,
    getdualvalue,
    getsolutiontime, 
    blockids

# DspModel placeholder
model = DspModel()

###############################################################################
# Override JuMP.Model
###############################################################################

# This is for the master problem.
function JuMP.Model(nblocks::Integer, master_has_subblocks::Bool = false)
    # construct model
    m = JuMP.Model()
    # free DspModel
    DspCInterface.freeModel(Dsp.model)
    # set block ids
    DspCInterface.setBlockIds(Dsp.model, nblocks, master_has_subblocks)
    # set extension
    m.ext[:DspBlocks] = BlockStructure(nothing, Dict{Int,JuMP.Model}(), Dict{Int,Float64}())
    # set solvehook
    JuMP.setsolvehook(m, Dsp.dsp_solve)
    # return model
    m
end

# This is for the subproblems.
function JuMP.Model(parent::JuMP.Model, id::Integer, weight::Float64)
    # construct model
    m = JuMP.Model()
    # set extension
    m.ext[:DspBlocks] = BlockStructure(nothing, Dict{Int,JuMP.Model}(), Dict{Int,Float64}())
    # set block
    m.ext[:DspBlocks].parent = parent
    setindex!(parent.ext[:DspBlocks].children, m, id)
    setindex!(parent.ext[:DspBlocks].weight, weight, id)
    # return model
    m
end

###############################################################################
# JuMP.solvehook
###############################################################################

# This function is hooked by JuMP (see block.jl)
function dsp_solve(m::JuMP.Model; suppress_warnings = false, options...)

    setoptions(options)

    # load JuMP model to Dsp
    DspCInterface.loadProblem(Dsp.model, m)

    # solve
    DspCInterface.solve(Dsp.model)

    # solution status
    statcode = DspCInterface.getStatus(Dsp.model)
    stat = parseStatusCode(statcode)
    
    if Dsp.model.solve_type != :DW
        # Extract solution from the solver
        Dsp.model.numRows = DspCInterface.getTotalNumRows(Dsp.model)
        Dsp.model.numCols = DspCInterface.getTotalNumCols(Dsp.model)
        m.objVal = NaN
        m.colVal = fill(NaN, Dsp.model.numCols)

        if stat != :Optimal
            suppress_warnings || warn("Not solved to optimality, status: $stat")
        end

        if !(stat == :Infeasible || stat == :Unbounded)
            getDspSolution(m)
        end
    end
    
    # Return the solve status
    stat
end

function setoptions(options)
    for (optname, optval) in options
        if optname == :param
            DspCInterface.readParamFile(Dsp.model, optval)
        elseif optname == :solve_type
            if optval in [:Dual, :Benders, :Extensive, :DW]
                Dsp.model.solve_type = optval
            else
                warn("solve_type $optval is not available.")
            end
        else
            warn("Options $optname is not available.")
        end
    end
end

###############################################################################
# Input/output files
###############################################################################

function optimize(;suppress_warnings = false, options...)

    setoptions(options)

    # solve
    DspCInterface.solve(Dsp.model)

    # solution status
    statcode = DspCInterface.getStatus(Dsp.model)
    stat = parseStatusCode(statcode)

    # Extract solution from the solver
    Dsp.model.numRows = DspCInterface.getTotalNumRows(Dsp.model)
    Dsp.model.numCols = DspCInterface.getTotalNumCols(Dsp.model)

    if stat != :Optimal
        suppress_warnings || warn("Not solved to optimality, status: $stat")
    end
#=
    if !(stat == :Infeasible || stat == :Unbounded)
        getDspSolution(nothing)
    end
=#
    # Return the solve status
    stat
end
JuMP.solve(;suppress_warnings = false, options...) = optimize(;suppress_warnings = false, options...)

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
# get primal value
getprimvalue() = Dsp.model.colVal
# get dual value
getdualvalue() = Dsp.model.rowVal
# get solution time
getsolutiontime() = DspCInterface.getWallTime(Dsp.model)
# get block ids
blockids() = Dsp.model.block_ids

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

function getDspSolution(m = nothing)
    Dsp.model.primVal = DspCInterface.getPrimalBound(Dsp.model)
    Dsp.model.dualVal = DspCInterface.getDualBound(Dsp.model)

    if Dsp.model.solve_type == :Dual
        Dsp.model.rowVal = DspCInterface.getDualSolution(Dsp.model)
    else
        Dsp.model.colVal = DspCInterface.getSolution(Dsp.model)
        if m != nothing
            # parse solution to each block
            n_start = 1
            n_end = m.numCols
            m.colVal = Dsp.model.colVal[n_start:n_end]
            n_start += m.numCols
            if haskey(m.ext, :DspBlocks) == true
                numBlockCols = DspCInterface.getNumBlockCols(Dsp.model, m)
                blocks = m.ext[:DspBlocks].children
                for i in 1:Dsp.model.nblocks
                    n_end += numBlockCols[i]
                    if haskey(blocks, i)
                        # @show b
                        # @show n_start
                        # @show n_end
                        blocks[i].colVal = Dsp.model.colVal[n_start:n_end]
                    end
                    n_start += numBlockCols[i]
                end
            end
        end
    end
    if m != nothing
        m.objVal = Dsp.model.primVal
        # maximization?
        if m.objSense == :Max
            m.objVal *= -1
            Dsp.model.primVal *= -1
            Dsp.model.dualVal *= -1
        end
    end
end

end # module
