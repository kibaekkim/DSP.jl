function readSmps(prob::DspModel, filename::AbstractString, dedicatedMaster::Bool)
    # Check pointer to TssModel
    check_problem()
    @dsp_ccall("readSmps", Void, (Ptr{Void}, Ptr{UInt8}), prob.p, convert(Vector{UInt8}, filename))
    nscen = getNumScenarios();
    proc_idx_set = 1:nscen;
    if isdefined(:MPI) == true
        proc_idx_set = getProcIdxSet(nscen, dedicatedMaster);
    end
    setProcIdxSet(prob, proc_idx_set);  
end
readSmps(prob::DspModel, filename::AbstractString) = readSmps(prob, filename, false);

function loadProblem(model::JuMP.Model, dedicatedMaster::Bool)
    check_problem()
    if haskey(model.ext, :DspBlocks)
        loadStochasticProblem(model, dedicatedMaster);
    else
        warn("No block is defined.")
        loadDeterministicProblem(model);
    end
end
loadProblem(model::JuMP.Model) = loadProblem(model, false);

function loadStochasticProblem(model::JuMP.Model, dedicatedMaster::Bool)
    
    # retrieve DspMathProgModel
    prob = model.internalModel

    # get DspBlocks
    blocks = model.ext[:DspBlocks]
    
    nscen  = convert(Cint, length(blocks.block))
    ncols1 = convert(Cint, model.numCols)
    nrows1 = convert(Cint, length(model.linconstr))
    ncols2 = 0
    nrows2 = 0
    proc_idx_set = collect(1:nscen);
    if isdefined(:MPI) == true && MPI.Initialized() == true
        proc_idx_set = getProcIdxSet(nscen, dedicatedMaster);
    end
    setProcIdxSet(proc_idx_set);
    for s in 1:length(proc_idx_set)
        ncols2 = convert(Cint, blocks.block[s].numCols)
        nrows2 = convert(Cint, length(blocks.block[s].linconstr))
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
        blk = blocks.block[s]
        probability = blocks.weight[s]
        # get model data
        start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(blk)
        @dsp_ccall("loadSecondStage", Void, 
            (Ptr{Void}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, 
                Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), 
            prob.p, proc_idx_set[s]-1, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)
    end
    
end

function loadDeterministicProblem(model::JuMP.Model)
    # retrieve DspMathProgModel
    prob = model.internalModel
    ncols = convert(Cint, model.numCols)
    nrows = convert(Cint, length(model.linconstr))
    start, index, value, clbd, cubd, ctype, obj, rlbd, rubd = getDataFormat(model)
    numels = length(index)
    @dsp_ccall("loadDeterministic", Void, 
        (Ptr{Void}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint, Cint, 
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            prob.p, start, index, value, numels, ncols, nrows, clbd, cubd, ctype, obj, rlbd, rubd)
end