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
    solver
    function DspModel()
        p = @dsp_ccall("createEnv", Ptr{Void}, ())
        solver = :DualDecomp
        prob = new(p)
        finalizer(prob, freeDSP)
        return prob
    end
end

function freeDSP(prob::DspModel)
    if prob.p == C_NULL
        return
    end
    @dsp_ccall("freeEnv", Void, (Ptr{Void},), prob.p)
    prob.p = C_NULL
    return
end

function check_problem(prob::DspModel)
    if prob.p == C_NULL
        error("Invalid DspModel")
    end
    return true
end

function readParamFile(prob::DspModel, param_file::AbstractString)
    @dsp_ccall("readParamFile", Void, (Ptr{Void}, Ptr{UInt8}), prob.p, param_file);
end

function prepConstrMatrix(m::JuMP.Model)
    if !haskey(m.ext, :DspBlocks)
        return JuMP.prepConstrMatrix(m)
    end

    blocks = m.ext[:DspBlocks]
    if stoch.parent == nothing
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
                if m.linconstr[nrow].terms.vars[id].m == stoch.parent
                    push!(cind, var.col)
                elseif m.linconstr[nrow].terms.vars[id].m == m
                    push!(cind, stoch.parent.numCols + var.col)
                end
                push!(value, aff.coeffs[id])
                splice!(aff.vars, id)
                splice!(aff.coeffs, id)
            end
        end
    end
    return sparse(rind, cind, value, length(m.linconstr), stoch.parent.numCols + m.numCols)
end


