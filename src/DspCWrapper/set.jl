function setSolverType(prob::DspModel, solver)
    solver_types = [:DualDecomp, :Benders, :ExtensiveForm]
    if solver in solver_types
        prob.solver = solver
    else
        warn("Solver type $solver is invalid.")
    end
end

function setProcIdxSet(prob::DspModel, scenarios::Array{Int,1})
    num = convert(Cint, length(scenarios));
    scenarios = scenarios - 1;
    setIntPtrParam(prob, "ARR_PROC_IDX", num, scenarios);
end

function setIntPtrParam(prob::DspModel, name::AbstractString, size::Integer, value::Array{Int,1})
    @dsp_ccall("setIntPtrParam", Void, (Ptr{Void}, Ptr{UInt8}, Cint, Ptr{Cint}), 
        prob.p, name, convert(Cint, size), convert(Vector{Cint}, value))
end