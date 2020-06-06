module Dsp

@enum Methods begin
    ExtensiveForm
    Benders
    Dual
    Legacy
end

mutable struct Model
    p::Ptr{Cvoid}

    # Number of blocks
    nblocks::Int

    # solve_type should be one of these:
    solve_type::Methods

    numRows::Int
    numCols::Int
    primVal::Float64
    dualVal::Float64
    colVal::Vector{Float64}
    rowVal::Vector{Float64}

    # MPI settings
    comm
    comm_size::Int
    comm_rank::Int

    # Array of block ids:
    # The size of array is not necessarily same as nblocks,
    # as block ids may be distributed to multiple processors.
    block_ids::Vector{Integer}

    is_stochastic::Bool
end

include("DspCInterface.jl")
using .DspCInterface

using Libdl
using StructJuMP
using SparseArrays
import MathOptInterface

const SJ = StructJuMP
const MOI = MathOptInterface

# Initialize Dsp.Model
dsp_model = Model(C_NULL, 0, Dual, 0, 0, NaN, NaN, [], [], nothing, 1, 0, [], false)

function __init__()
    try
        # path_to_lib = ENV["JULIA_DSP_LIBRARY_PATH"]
        if Sys.islinux()
            Libdl.dlopen("libDsp.so", Libdl.RTLD_GLOBAL)
        elseif Sys.isapple()
            Libdl.dlopen("libDsp.dylib", Libdl.RTLD_GLOBAL)
        end
        dsp_model.p = DspCInterface.createEnv()
        finalizer(DspCInterface.freeEnv, dsp_model)
    catch
        @warn("Could not load DSP shared library. Make sure it is in your library path.")
        rethrow()
    end
end

freeSolver(dsp::Model) = DspCInterface.freeSolver(dsp)

function freeModel(dsp::Model)
    @assert(dsp.p != C_NULL)
    DspCInterface.freeModel(dsp)
    dsp.nblocks = 0
    dsp.solve_type = Dual
    dsp.numRows = 0
    dsp.numCols = 0
    dsp.primVal = NaN
    dsp.dualVal = NaN
    dsp.colVal = Vector{Float64}()
    dsp.rowVal = Vector{Float64}()
    dsp.block_ids = Vector{Integer}()
    dsp.is_stochastic = false
end

function freeEnv(dsp::Model)
    @assert(dsp.p != C_NULL)
    freeModel(dsp)
    DspCInterface.freeEnv(dsp)
    dsp.p = C_NULL
    dsp.comm = nothing
    dsp.comm_size = 1
    dsp.comm_rank = 0
end

function optimize!(m::SJ.StructuredModel; options...)
    @assert(dsp_model.p != C_NULL)

    # remove existing model
    DspCInterface.freeModel(dsp_model)

    setoptions(dsp_model, options)
    loadProblem(dsp_model, m)
    solve(dsp_model)

    # solution status
    statcode = DspCInterface.getStatus(dsp_model)
end

function readSmps(dsp::Model, filename::AbstractString, master_has_subblocks::Bool = false)
    @assert(dsp.p != C_NULL)
    DspCInterface.readSmps(dsp, filename)
    setBlockIds(dsp, DspCInterface.getNumScenarios(dsp), master_has_subblocks)
end

function get_model_data(m::SJ.StructuredModel)

    # Get a column-wise sparse matrix
    start, index, value, rlbd, rubd = get_constraint_matrix(m)

    # column information
    clbd = Vector{Float64}(undef, num_variables(m))
    cubd = Vector{Float64}(undef, num_variables(m))
    ctype = ""
    cname = Vector{String}(undef, num_variables(m))
    for i in 1:num_variables(m)
        vref = SJ.StructuredVariableRef(m, i)
        v = m.variables[vref.idx]
        if v.info.integer
            ctype = ctype * "I"
        elseif v.info.binary
            ctype = ctype * "B"
        else
            ctype = ctype * "C"
        end
        clbd[vref.idx] = v.info.has_lb ? v.info.lower_bound : -Inf
        cubd[vref.idx] = v.info.has_ub ? v.info.upper_bound : Inf
        cname[vref.idx] = m.varnames[vref.idx]
    end

    # objective coefficients
    obj = zeros(num_variables(m))
    if !(objective_function_type(m) <: Real)
        for (v,coef) in objective_function(m).terms
            obj[v.idx] = coef
        end
    end

    if objective_sense(m) == MOI.MAX_SENSE
        obj .*= -1
    end

    return start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname
end

function get_constraint_matrix(m::SJ.StructuredModel)

    is_parent = SJ.getparent(m) == nothing ? true : false

    num_rows = 0 # need to count
    num_cols = num_variables(m)
    if !is_parent
        num_cols += num_variables(SJ.getparent(m))
    end

    # count the number of nonzero elements
    nnz = 0
    for (i,cons) in m.constraints
        nnz += length(cons.func.terms)
        num_rows += 1
    end

    rind = Vector{Int}(undef, nnz)
    cind = Vector{Int}(undef, nnz)
    value = Vector{Float64}(undef, nnz)
    rlbd = Vector{Float64}(undef, num_rows)
    rubd = Vector{Float64}(undef, num_rows)

    pos = 1
    for (i,cons) in m.constraints
        for (var,coef) in cons.func.terms
            rind[pos] = i
            if is_parent
                cind[pos] = var.idx
            elseif JuMP.owner_model(var) == SJ.getparent(m)
                cind[pos] = var.idx
            else
                cind[pos] = var.idx + num_variables(SJ.getparent(m))
            end
            value[pos] = coef
            pos += 1
        end

        # row bounds
        rlbd[i], rubd[i] = row_bounds_from_moi(cons.set)
    end
    @assert(pos-1==nnz)

    # NOTE: DSP takes CSR (row-wise sparse matrix) format.
    # So I simply switch the entries for columns and rows.
    mat = sparse(cind, rind, value, num_cols, num_rows)
    dropzeros!(mat)

    # sparse description
    start = convert(Vector{Cint}, mat.colptr .- 1)
    index = convert(Vector{Cint}, mat.rowval .- 1)
    value = mat.nzval

    return start, index, value, rlbd, rubd
end

row_bounds_from_moi(set::MOI.LessThan) = (-Inf, set.upper)
row_bounds_from_moi(set::MOI.GreaterThan) = (set.lower, Inf)
row_bounds_from_moi(set::MOI.EqualTo) = (set.value, set.value)
row_bounds_from_moi(set::MOI.Interval) = @error("Interval row bounds are not supported.")

function setoptions(dsp::Model, options)
    for (optname, optval) in options
        if optname == :param
            DspCInterface.readParamFile(dsp, optval)
        elseif optname == :is_stochastic
            dsp.is_stochastic = optval
        elseif optname == :solve_type
            if optval in instances(Methods)
                dsp.solve_type = optval
            else
                @warn("solve_type $optval is not available.")
            end
        else
            @warn("Options $optname is not available.")
        end
    end
end

function loadProblem(dsp::Model, model::SJ.StructuredModel)
    DspCInterface.freeModel(dsp)
    if dsp.is_stochastic
        loadStochasticProblem(dsp, model)
    else
        loadStructuredProblem(dsp, model)
    end
end

function loadStochasticProblem(dsp::Model, model::SJ.StructuredModel)

    nscen = SJ.num_scenarios(model)
    ncols1 = length(model.variables)
    nrows1 = length(model.constraints)
    ncols2 = 0
    nrows2 = 0
    for subm in values(SJ.getchildren(model))
        ncols2 = length(subm.variables)
        nrows2 = length(subm.constraints)
        break
    end

    # set scenario indices for each MPI processor
    if dsp.comm_size > 1
        ncols2 = MPI.allreduce([ncols2], MPI.MAX, dsp.comm)[1]
        nrows2 = MPI.allreduce([nrows2], MPI.MAX, dsp.comm)[1]
    end

    DspCInterface.setNumberOfScenarios(dsp, nscen)
    DspCInterface.setDimensions(dsp, ncols1, nrows1, ncols2, nrows2)

    # get problem data
    start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = get_model_data(model)
    DspCInterface.loadFirstStage(dsp, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

    for (id, blk) in SJ.getchildren(model)
        probability = SJ.getprobability(model)[id]
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = get_model_data(blk)
        DspCInterface.loadSecondStage(dsp, id-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
    end
end

function loadStructuredProblem(dsp::Model, model::SJ.StructuredModel)

    ncols1 = length(model.variables)
    nrows1 = length(model.constraints)

    # TODO: do something for MPI
    
    # load master
    start, index, value, rlbd, rubd, obj1, clbd1, cubd1, ctype1, cname = get_model_data(model)
    DspCInterface.loadBlockProblem(dsp, 0, ncols1, nrows1, start[end],
        start, index, value, clbd1, cubd1, ctype1, obj1, rlbd, rubd)

    # going over blocks
    for (id, blk) in SJ.getchildren(model)
        probability = SJ.getprobability(model)[id]
        ncols2 = length(blk.variables)
        nrows2 = length(blk.constraints)
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = get_model_data(blk)
        DspCInterface.loadBlockProblem(dsp, id, ncols1 + ncols2, nrows2, start[end], 
            start, index, value, [clbd1; clbd], [cubd1; cubd], [ctype1; ctype], [obj1; obj], rlbd, rubd)
    end

    # Finalize loading blocks
    DspCInterface.updateBlocks(dsp)
end

function solve(dsp::Model)
    if dsp.comm_size == 1
        if dsp.solve_type == Legacy
            if dsp.is_stochastic
                DspCInterface.solveDd(dsp);
            else
                @warn("This method is available for stochastic programs only.")
            end
        elseif dsp.solve_type == Benders
            if dsp.is_stochastic
                DspCInterface.solveBd(dsp);
            else
                @warn("This method is available for stochastic programs only.")
            end
        elseif dsp.solve_type == ExtensiveForm
            DspCInterface.solveDe(dsp);
        elseif dsp.solve_type == Dual
            DspCInterface.solveDw(dsp);
        end
    elseif dsp.comm_size > 1
        if dsp.solve_type == Legacy
            if dsp.is_stochastic
                DspCInterface.solveDdMpi(dsp, dsp.comm);
            else
                @warn("This method is available for stochastic programs only.")
            end
        elseif dsp.solve_type == Benders
            if dsp.is_stochastic
                DspCInterface.solveBdMpi(dsp, dsp.comm);
            else
                @warn("This method is available for stochastic programs only.")
            end
        elseif dsp.solve_type == Dual
            DspCInterface.solveDwMpi(dsp, dsp.comm);
        elseif dsp.solve_type == ExtensiveForm
            DspCInterface.solveDe(dsp);
        end
    end
end

function setBlockIds(dsp::Model, nblocks::Int, master_has_subblocks::Bool = false)
    # set number of blocks
    dsp.nblocks = nblocks
    # set MPI settings
    if @isdefined(MPI) && MPI.Initialized()
        dsp.comm = MPI.COMM_WORLD
        dsp.comm_size = MPI.Comm_size(dsp.comm)
        dsp.comm_rank = MPI.Comm_rank(dsp.comm)
    end
    #@show dsp.nblocks
    #@show dsp.comm
    #@show dsp.comm_size
    #@show dsp.comm_rank
    # get block ids with MPI settings
    dsp.block_ids = getBlockIds(dsp, master_has_subblocks)
    #@show dsp.block_ids
    # send the block ids to Dsp
    DspCInterface.setIntPtrParam(dsp, "ARR_PROC_IDX", length(dsp.block_ids), dsp.block_ids .- 1)
end

function getBlockIds(dsp::Model, master_has_subblocks::Bool = false)
    # processor info
    mysize = dsp.comm_size
    myrank = dsp.comm_rank
    # empty block ids
    proc_idx_set = Int[]
    # DSP is further parallelized with mysize > dsp.nblocks.
    modrank = myrank % dsp.nblocks
    # If we have more than one processor,
    # do not assign a sub-block to the master.
    if master_has_subblocks
        # assign sub-blocks in round-robin fashion
        for s = modrank:mysize:(dsp.nblocks-1)
            push!(proc_idx_set, s+1)
        end
    else
        if mysize > 1
            if myrank == 0
                return proc_idx_set
            end
            # exclude master
            mysize -= 1;
            modrank = (myrank-1) % dsp.nblocks
        end
        # assign sub-blocks in round-robin fashion
        for s = modrank:mysize:(dsp.nblocks-1)
            push!(proc_idx_set, s+1)
        end
    end
    # return assigned block ids
    return proc_idx_set
end

function getNumBlockCols(dsp::Model, model::SJ.StructuredModel)
    # subblocks
    blocks = SJ.getchildren(model)
    # get number of block columns
    numBlockCols = Dict{Int,Int}()
    if dsp.comm_size > 1
        num_proc_blocks = convert(Vector{Cint}, MPI.Allgather(length(blocks), dsp.comm))
        #@show num_proc_blocks
        #@show collect(keys(blocks))
        block_ids = MPI.Allgatherv(collect(keys(blocks)), num_proc_blocks, dsp.comm)
        #@show block_ids
        ncols_to_send = Int[blocks[i].numCols for i in keys(blocks)]
        #@show ncols_to_send
        ncols = MPI.Allgatherv(ncols_to_send, num_proc_blocks, dsp.comm)
        #@show ncols
        for i in 1:dsp.nblocks
            setindex!(numBlockCols, ncols[i], block_ids[i])
        end
    else
        for b in blocks
            setindex!(numBlockCols, b.second.numCols, b.first)
        end
    end
    return numBlockCols
end

#=

import JuMP
export
    readSmps, 
    writeMps,
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

mutable struct BlockStructure
    parent
    children::Dict{Int,JuMP.Model}
    weight::Dict{Int,Float64}
end

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
    m.optimize_hook = Dsp.dsp_solve
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
    # parse options
    setoptions(options)

    # load JuMP model to Dsp
    DspCInterface.loadProblem(Dsp.model, m)

    # solve
    DspCInterface.solve(Dsp.model)

    # solution status
    statcode = DspCInterface.getStatus(Dsp.model)
    stat = parseStatusCode(statcode)
    
    # Extract solution from the solver
    Dsp.model.numRows = DspCInterface.getTotalNumRows(Dsp.model)
    Dsp.model.numCols = DspCInterface.getTotalNumCols(Dsp.model)
    m.objVal = NaN
    m.colVal = fill(NaN, Dsp.model.numCols)

    if stat != :Optimal
        suppress_warnings || @warn("Not solved to optimality, status: $stat")
    end

    if !(stat == :Infeasible || stat == :Unbounded)
        getDspSolution(m)
    end
    
    # Return the solve status
    stat
end

###############################################################################
# Input/output files
###############################################################################

function optimize(;suppress_warnings = false, options...)
    # parse options
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
        suppress_warnings || @warn("Not solved to optimality, status: $stat")
    end

    if Dsp.model.solve_type != :BB
        if !(stat == :Infeasible || stat == :Unbounded)
            getDspSolution()
        end
    end
    # Return the solve status
    stat
end
JuMP.solve(;suppress_warnings = false, options...) = optimize(;suppress_warnings = false, options...)

# Read model from SMPS files
function readSmps(filename::AbstractString, master_has_subblocks::Bool = false)
    # free DspModel
    DspCInterface.freeModel(Dsp.model)
    # read Smps file
    DspCInterface.readSmps(Dsp.model, filename, master_has_subblocks)
end

# Write model to MPS file
function writeMps(filename::AbstractString)
    DspCInterface.@dsp_ccall("writeMps", Cvoid, (Ptr{Cvoid}, Ptr{UInt8}), Dsp.model.p, filename)
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
        @warn("No block was created in the model.")
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
        @warn("Unknown status: $statcode")
    end

    stat
end

function getDspSolution(m = nothing)
    Dsp.model.primVal = DspCInterface.getPrimalBound(Dsp.model)
    Dsp.model.dualVal = DspCInterface.getDualBound(Dsp.model)

    if Dsp.model.solve_type == :Dual
        Dsp.model.rowVal = DspCInterface.getDualSolution(Dsp.model)
	end

	if Dsp.model.solve_type != :BB || Dsp.model.primVal < 1.0e+20
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
    
            m.objVal = Dsp.model.primVal
            # maximization?
            if m.objSense == :Max
                m.objVal *= -1
                Dsp.model.primVal *= -1
                Dsp.model.dualVal *= -1
            end
        end
	end
end
=#
end # module
