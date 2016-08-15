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

    function DspModel()
        p = @dsp_ccall("createEnv", Ptr{Void}, ())
        solve_type = :Dual
        numRows = 0
        numCols = 0
        primVal = NaN
        dualVal = NaN
        colVal = Vector{Float64}()
        rowVal = Vector{Float64}()
        prob = new(p, solve_type, numRows, numCols, primVal, dualVal, colVal, rowVal)
        finalizer(prob, freeDSP)
        return prob
    end
end

function freeDSP(prob::DspModel)
    if prob.p == C_NULL
        return
    else
        @dsp_ccall("freeEnv", Void, (Ptr{Void},), prob.p)
        prob.p = C_NULL
    end
    prob.solve_type = nothing
    prob.numRows = 0
    prob.numCols = 0
    prob.primVal = NaN
    prob.dualVal = NaN
    prob.colVal = Vector{Float64}()
    prob.rowVal = Vector{Float64}()
    return
end

function freeModel(prob::DspModel)
    check_problem(prob)
    @dsp_ccall("freeTssModel", Void, (Ptr{Void},), prob.p)
end

function check_problem(prob::DspModel)
    if prob.p == C_NULL
        error("Invalid DspModel")
    end
    return true
end

function readParamFile(prob::DspModel, param_file::AbstractString)
    check_problem(prob)
    @dsp_ccall("readParamFile", Void, (Ptr{Void}, Ptr{UInt8}), prob.p, param_file);
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
# Load problems
###############################################################################

function readSmps(prob::DspModel, filename::AbstractString, dedicatedMaster::Bool)
    # Check pointer to TssModel
    check_problem(prob)
    @dsp_ccall("readSmps", Void, (Ptr{Void}, Ptr{UInt8}), prob.p, convert(Vector{UInt8}, filename))
    nscen = getNumScenarios(prob);
    proc_idx_set = collect(1:nscen);
    if isdefined(:MPI) == true
        proc_idx_set = getProcIdxSet(nscen, dedicatedMaster);
    end
    setProcIdxSet(prob, proc_idx_set);  
end
readSmps(prob::DspModel, filename::AbstractString) = readSmps(prob, filename, false);

function loadProblem(prob::DspModel, model::JuMP.Model, dedicatedMaster::Bool)
    check_problem(prob)
    if haskey(model.ext, :DspBlocks)
        loadStochasticProblem(prob, model, dedicatedMaster);
    else
        warn("No block is defined.")
        loadDeterministicProblem(prob, model);
    end
end
loadProblem(prob::DspModel, model::JuMP.Model) = loadProblem(prob, model, false);

function loadStochasticProblem(prob::DspModel, model::JuMP.Model, dedicatedMaster::Bool)
    
    # get DspBlocks
    blocks = model.ext[:DspBlocks]
    
    nscen  = convert(Cint, length(blocks.children))
    ncols1 = convert(Cint, model.numCols)
    nrows1 = convert(Cint, length(model.linconstr))
    ncols2 = 0
    nrows2 = 0
    proc_idx_set = collect(1:nscen);
    if isdefined(:MPI) == true && MPI.Initialized() == true
        proc_idx_set = getProcIdxSet(nscen, dedicatedMaster);
    end
    setProcIdxSet(prob, proc_idx_set);
    for s in 1:length(proc_idx_set)
        ncols2 = convert(Cint, blocks.children[s].numCols)
        nrows2 = convert(Cint, length(blocks.children[s].linconstr))
        break;
    end
    
    @dsp_ccall("setNumberOfScenarios", Void, (Ptr{Void}, Cint), prob.p, nscen)
    @dsp_ccall("setDimensions", Void, 
        (Ptr{Void}, Cint, Cint, Cint, Cint), 
        prob.p, ncols1, nrows1, ncols2, nrows2)
    
    # get problem data
    start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(model)
    
    @dsp_ccall("loadFirstStage", Void, 
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, 
            Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            prob.p, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
    
    for s in 1:length(proc_idx_set)
        # get model
        blk = blocks.children[s]
        probability = blocks.weight[s]
        # get model data
        start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(blk)
        @dsp_ccall("loadSecondStage", Void, 
            (Ptr{Void}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, 
                Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
            prob.p, proc_idx_set[s]-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
    end
    
end

function loadDeterministicProblem(prob::DspModel, model::JuMP.Model)
    ncols = convert(Cint, model.numCols)
    nrows = convert(Cint, length(model.linconstr))
    start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(model)
    numels = length(index)
    @dsp_ccall("loadDeterministic", Void, 
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint, Cint, 
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            prob.p, start, index, value, numels, ncols, nrows, clbd, cubd, ctype, obj, rlbd, rubd)
end


###############################################################################
# Get functions
###############################################################################

function solveDe(prob::DspModel)
    check_problem(prob)
    return @dsp_ccall("solveDe", Void, (Ptr{Void},), prob.p)
end

function solveBd(prob::DspModel)
    check_problem(prob)
    return @dsp_ccall("solveBd", Void, (Ptr{Void},), prob.p)
end

function solveDd(prob::DspModel, comm)
    check_problem(prob)
    return @dsp_ccall("solveDd", Void, (Ptr{Void}, Cint), prob.p, convert(Cint, comm.val))
end

function solveBdMpi(prob::DspModel, comm)
    check_problem(prob)
    return @dsp_ccall("solveBdMpi", Void, (Ptr{Void}, Cint, Cint), prob.p, 1, convert(Cint, comm.val))
end

###############################################################################
# Get functions
###############################################################################

function getProcIdxSet(numScens::Integer, dedicatedMaster::Bool)
    mysize = 1;
    myrank = 0;
    if isdefined(:MPI) == true && MPI.Initialized() == true
        comm = MPI.COMM_WORLD
        mysize = MPI.Comm_size(comm)
        myrank = MPI.Comm_rank(comm)
    end
    # Round-and-Robin
    proc_idx_set = Int[];
    # DSP is further parallelized with mysize > numScens.
    modrank = myrank % numScens;
    
    if dedicatedMaster == true
        if myrank == 0
            return proc_idx_set;
        end
        # exclude master
        mysize -= 1;
        modrank = (myrank-1) % numScens;
    end
    
    for s = modrank:mysize:(numScens-1)
        push!(proc_idx_set, s+1);
    end
    return proc_idx_set;
end
getProcIdxSet(numScens::Integer) = getProcIdxSet(numScens,false);

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
                   (:getSolutionStatus, Cint), 
                   (:getNumIterations, Cint), 
                   (:getNumNodes, Cint), 
                   (:getSolutionTime, Cdouble), 
                   (:getPrimalBound, Cdouble), 
                   (:getDualBound, Cdouble)]
    strfunc = string(func)
    @eval begin
        function $func(prob::DspModel)
            check_problem(prob)
            return @dsp_ccall($strfunc, $rtn, (Ptr{Void},), prob.p)
        end
    end
end

function getObjCoef(prob::DspModel)
    check_problem(prob)
    num = getTotalNumCols()
    obj = Array(Cdouble, num)
    @dsp_ccall("getObjCoef", Void, (Ptr{Void}, Ptr{Cdouble}), prob.p, obj)
    return obj
end

function getSolution(prob::DspModel, num::Integer)
    sol = Array(Cdouble, num)
    @dsp_ccall("getSolution", Void, (Ptr{Void}, Cint, Ptr{Cdouble}), prob.p, num, sol)
    return sol
end
getSolution(prob::DspModel) = getSolution(prob, getTotalNumCols(prob))

###############################################################################
# Set functions
###############################################################################

function setSolverType(prob::DspModel, solver)
    check_problem(prob)
    solver_types = [:DualDecomp, :Benders, :ExtensiveForm]
    if solver in solver_types
        prob.solver = solver
    else
        warn("Solver type $solver is invalid.")
    end
end

function setProcIdxSet(prob::DspModel, scenarios::Array{Int,1})
    check_problem(prob)
    num = convert(Cint, length(scenarios));
    scenarios = scenarios - 1;
    @dsp_ccall("setIntPtrParam", Void, (Ptr{Void}, Ptr{UInt8}, Cint, Ptr{Cint}), 
        prob.p, "ARR_PROC_IDX", convert(Cint, num), convert(Vector{Cint}, scenarios))
end

end # end of module