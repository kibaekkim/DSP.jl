module DspCInterface

import ..Dsp
import Compat: String, unsafe_wrap
import JuMP

Pkg.installed("MPI") == nothing || using MPI

export DspModel

###############################################################################
# Help functions
###############################################################################

macro dsp_ccall(func, args...)
    @unix_only return quote
        ccall(($func, "libDsp"), $(args...))
    end
    @windows_only return quote
        ccall(($func, "libDsp"), stdcall, $(args...))
    end
end

type DspModel
    p::Ptr{Void}

    # Number of blocks
    nblocks::Int

    # solve_type should be one of these:
    # :Dual
    # :Benders
    # :Extensive
    solve_type

    numRows::Int
    numCols::Int
    primVal
    dualVal
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

    function DspModel()
        # assign Dsp pointer
        p = @dsp_ccall("createEnv", Ptr{Void}, ())
        # initialize variables
        nblocks = 0
        solve_type = :Dual
        numRows = 0
        numCols = 0
        primVal = NaN
        dualVal = NaN
        colVal = Vector{Float64}()
        rowVal = Vector{Float64}()
        comm = nothing
        comm_size = 1
        comm_rank = 0
        block_ids = Vector{Integer}()
        # create DspModel
        dsp = new(p, nblocks, solve_type, numRows, numCols, primVal, dualVal, colVal, rowVal, comm, comm_size, comm_rank, block_ids)
        # with finalizer
        finalizer(dsp, freeDSP)
        # return DspModel
        return dsp
    end
end

function freeDSP(dsp::DspModel)
    if dsp.p == C_NULL
        return
    else
        @dsp_ccall("freeEnv", Void, (Ptr{Void},), dsp.p)
        dsp.p = C_NULL
    end
    dsp.nblocks = 0
    dsp.solve_type = nothing
    dsp.numRows = 0
    dsp.numCols = 0
    dsp.primVal = NaN
    dsp.dualVal = NaN
    dsp.colVal = Vector{Float64}()
    dsp.rowVal = Vector{Float64}()
    return
end

function freeModel(dsp::DspModel)
    check_problem(dsp)
    @dsp_ccall("freeModel", Void, (Ptr{Void},), dsp.p)
    dsp.nblocks = 0
    dsp.numRows = 0
    dsp.numCols = 0
    dsp.primVal = NaN
    dsp.dualVal = NaN
    dsp.colVal = Vector{Float64}()
    dsp.rowVal = Vector{Float64}()
end

function check_problem(dsp::DspModel)
    if dsp.p == C_NULL
        error("Invalid DspModel")
    end
end

function readParamFile(dsp::DspModel, param_file::AbstractString)
    check_problem(dsp)
    @dsp_ccall("readParamFile", Void, (Ptr{Void}, Ptr{UInt8}), dsp.p, param_file);
end

function prepConstrMatrix(m::JuMP.Model)
    if !haskey(m.ext, :DspBlocks)
        return JuMP.prepConstrMatrix(m)
    end

    blocks = m.ext[:DspBlocks]
    if blocks.parent == nothing
        return JuMP.prepConstrMatrix(m)
    else
        rind = Int[]
        cind = Int[]
        value = Float64[]
        linconstr = deepcopy(m.linconstr)
        for (nrow,con) in enumerate(linconstr)
            aff = con.terms
            for (var,id) in zip(reverse(aff.vars), length(aff.vars):-1:1)
                push!(rind, nrow)
                if m.linconstr[nrow].terms.vars[id].m == blocks.parent
                    push!(cind, var.col)
                elseif m.linconstr[nrow].terms.vars[id].m == m
                    push!(cind, blocks.parent.numCols + var.col)
                end
                push!(value, aff.coeffs[id])
                splice!(aff.vars, id)
                splice!(aff.coeffs, id)
            end
        end
    end
    return sparse(rind, cind, value, length(m.linconstr), blocks.parent.numCols + m.numCols)
end

###############################################################################
# Block IDs
###############################################################################

function setBlockIds(dsp::DspModel, nblocks::Integer)
    check_problem(dsp)
    # set number of blocks
    dsp.nblocks = nblocks
    # set MPI settings
    if isdefined(:MPI) && MPI.Initialized()
        dsp.comm = MPI.COMM_WORLD
        dsp.comm_size = MPI.Comm_size(dsp.comm)
        dsp.comm_rank = MPI.Comm_rank(dsp.comm)
    end
    #@show dsp.nblocks
    #@show dsp.comm
    #@show dsp.comm_size
    #@show dsp.comm_rank
    # get block ids with MPI settings
    dsp.block_ids = getBlockIds(dsp)
    #@show dsp.block_ids
    # send the block ids to Dsp
    @dsp_ccall("setIntPtrParam", Void, (Ptr{Void}, Ptr{UInt8}, Cint, Ptr{Cint}), 
        dsp.p, "ARR_PROC_IDX", convert(Cint, length(dsp.block_ids)), convert(Vector{Cint}, dsp.block_ids - 1))
end

function getBlockIds(dsp::DspModel)
    check_problem(dsp)
    # processor info
    mysize = dsp.comm_size
    myrank = dsp.comm_rank
    # empty block ids
    proc_idx_set = Int[]
    # DSP is further parallelized with mysize > dsp.nblocks.
    modrank = myrank % dsp.nblocks
    # If we have more than one processor, 
    # do not assign a sub-block to the master.
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
    # return assigned block ids
    return proc_idx_set
end

function getNumBlockCols(dsp::DspModel, m::JuMP.Model)
    check_problem(dsp)
    # subblocks
    blocks = m.ext[:DspBlocks].children
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

###############################################################################
# Load problems
###############################################################################

function readSmps(dsp::DspModel, filename::AbstractString)
    # Check pointer to TssModel
    check_problem(dsp)
    # read smps files
    @dsp_ccall("readSmps", Void, (Ptr{Void}, Ptr{UInt8}), dsp.p, convert(Vector{UInt8}, filename))
    # set block Ids
    setBlockIds(dsp, getNumScenarios(dsp))
end

function loadProblem(dsp::DspModel, model::JuMP.Model, dedicatedMaster::Bool)
    check_problem(dsp)
    if haskey(model.ext, :DspBlocks)
        loadStochasticProblem(dsp, model, dedicatedMaster)
    else
        warn("No block is defined.")
        loadDeterministicProblem(dsp, model)
    end
end
loadProblem(dsp::DspModel, model::JuMP.Model) = loadProblem(dsp, model, true);

function loadStochasticProblem(dsp::DspModel, model::JuMP.Model, dedicatedMaster::Bool)
    
    # get DspBlocks
    blocks = model.ext[:DspBlocks]
    
    nscen  = dsp.nblocks
    ncols1 = model.numCols
    nrows1 = length(model.linconstr)
    ncols2 = 0
    nrows2 = 0
    for s in values(blocks.children)
        ncols2 = s.numCols
        nrows2 = length(s.linconstr)
        break
    end
    
    # set scenario indices for each MPI processor
    if dsp.comm_size > 1
        ncols2 = MPI.allreduce([ncols2], MPI.MAX, dsp.comm)[1]
        nrows2 = MPI.allreduce([nrows2], MPI.MAX, dsp.comm)[1]
    end
    
    @dsp_ccall("setNumberOfScenarios", Void, (Ptr{Void}, Cint), dsp.p, convert(Cint, nscen))
    @dsp_ccall("setDimensions", Void, 
        (Ptr{Void}, Cint, Cint, Cint, Cint), 
        dsp.p, convert(Cint, ncols1), convert(Cint, nrows1), convert(Cint, ncols2), convert(Cint, nrows2))
    
    # get problem data
    start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(model)
    
    @dsp_ccall("loadFirstStage", Void, 
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, 
            Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            dsp.p, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

    for id in dsp.block_ids
        # model and probability
        blk = blocks.children[id]
        probability = blocks.weight[id]
        # get model data
        start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(blk)
        @dsp_ccall("loadSecondStage", Void, 
            (Ptr{Void}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, 
                Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
            dsp.p, id-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
    end
    
end

function loadDeterministicProblem(dsp::DspModel, model::JuMP.Model)
    ncols = convert(Cint, model.numCols)
    nrows = convert(Cint, length(model.linconstr))
    start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(model)
    numels = length(index)
    @dsp_ccall("loadDeterministic", Void, 
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint, Cint, 
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            dsp.p, start, index, value, numels, ncols, nrows, clbd, cubd, ctype, obj, rlbd, rubd)
end


###############################################################################
# Get functions
###############################################################################

for func in [:freeSolver, 
             :solveDe, 
             :solveBd, 
             :solveDd]
    strfunc = string(func)
    @eval begin
        function $func(dsp::DspModel)
            return @dsp_ccall($strfunc, Void, (Ptr{Void},), dsp.p)
        end
    end
end

for func in [:solveBdMpi, :solveDdMpi]
    strfunc = string(func)
    @eval begin
        function $func(dsp::DspModel, comm)
            return @dsp_ccall($strfunc, Void, (Ptr{Void}, Cint), dsp.p, convert(Cint, comm.val))
        end
    end
end

function solve(dsp::DspModel)
    check_problem(dsp)
    if dsp.comm_size == 1
        if dsp.solve_type == :Dual
            solveDd(dsp);
        elseif dsp.solve_type == :Benders
            solveBd(dsp);
        elseif dsp.solve_type == :Extensive
            solveDe(dsp);
        end
    elseif dsp.comm_size > 1
        if dsp.solve_type == :Dual
            solveDdMpi(dsp, dsp.comm);
        elseif dsp.solve_type == :Benders
            solveBdMpi(dsp, dsp.comm);
        elseif dsp.solve_type == :Extensive
            solveDe(dsp);
        end
    end
end

###############################################################################
# Get functions
###############################################################################

function getDataFormat(model::JuMP.Model)
    # Get a column-wise sparse matrix
    mat = prepConstrMatrix(model)
    
    # Tranpose; now I have row-wise sparse matrix
    mat = mat'
    
    # sparse description
    start = convert(Vector{Cint}, mat.colptr - 1)
    index = convert(Vector{Cint}, mat.rowval - 1)
    value = mat.nzval
    
    # column type
    ctype = ""
    for i = 1:length(model.colCat)
        if model.colCat[i] == :Int
            ctype = ctype * "I";
        elseif model.colCat[i] == :Bin
            ctype = ctype * "B";
        else
            ctype = ctype * "C";
        end
    end
    ctype = convert(Vector{UInt8}, ctype)
    
    # objective coefficients
    obj, rlbd, rubd = JuMP.prepProblemBounds(model)

    # set objective sense
    if model.objSense == :Max
        obj *= -1
    end
    
    return start, index, value, model.colLower, model.colUpper, ctype, obj, rlbd, rubd
end

for (func,rtn) in [(:getNumScenarios, Cint), 
                   (:getTotalNumCols, Cint), 
                   (:getStatus, Cint), 
                   (:getNumIterations, Cint), 
                   (:getNumNodes, Cint), 
                   (:getWallTime, Cdouble), 
                   (:getPrimalBound, Cdouble), 
                   (:getDualBound, Cdouble)]
    strfunc = string(func)
    @eval begin
        function $func(dsp::DspModel)
            check_problem(dsp)
            return @dsp_ccall($strfunc, $rtn, (Ptr{Void},), dsp.p)
        end
    end
end
getSolutionStatus(dsp::DspModel) = getStatus(dsp)

function getObjCoef(dsp::DspModel)
    check_problem(dsp)
    num = getTotalNumCols()
    obj = Array(Cdouble, num)
    @dsp_ccall("getObjCoef", Void, (Ptr{Void}, Ptr{Cdouble}), dsp.p, obj)
    return obj
end

function getSolution(dsp::DspModel, num::Integer)
    sol = Array(Cdouble, num)
    @dsp_ccall("getPrimalSolution", Void, (Ptr{Void}, Cint, Ptr{Cdouble}), dsp.p, num, sol)
    return sol
end
getSolution(dsp::DspModel) = getSolution(dsp, getTotalNumCols(dsp))

function getDualSolution(dsp::DspModel, num::Integer)
    sol = Array(Cdouble, num)
    @dsp_ccall("getDualSolution", Void, (Ptr{Void}, Cint, Ptr{Cdouble}), dsp.p, num, sol)
    return sol
end
getDualSolution(dsp::DspModel) = getDualSolution(dsp, getNumCouplingRows(dsp))

###############################################################################
# Set functions
###############################################################################

function setSolverType(dsp::DspModel, solver)
    check_problem(dsp)
    solver_types = [:DualDecomp, :Benders, :ExtensiveForm]
    if solver in solver_types
        dsp.solver = solver
    else
        warn("Solver type $solver is invalid.")
    end
end

end # end of module
