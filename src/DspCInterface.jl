module DspCInterface

using Pkg
using SparseArrays
import ..Dsp
using StructJuMP

const SJ = StructJuMP

# Use MPI only if installed
deps = Pkg.dependencies()
for (uuid, dep) in deps
    dep.is_direct_dep || continue
    dep.version === nothing && continue
    if dep.name == "MPI"
        using MPI
    end
end

###############################################################################
# Help functions
###############################################################################

macro dsp_ccall(func, args...)
    @static if !Sys.iswindows()
        return esc(quote
            ccall(($func, "libDsp"), $(args...))
        end)
    end
    @static if Sys.iswindows()
        return esc(quote
            ccall(($func, "libDsp"), stdcall, $(args...))
        end)
    end
end

createEnv() = @dsp_ccall("createEnv", Ptr{Cvoid}, ())
freeEnv(dsp::Dsp.Model) = @dsp_ccall("freeEnv", Cvoid, (Ptr{Cvoid},), dsp.p)
freeModel(dsp::Dsp.Model) = @dsp_ccall("freeModel", Cvoid, (Ptr{Cvoid},), dsp.p)

readParamFile(dsp::Dsp.Model, param_file::AbstractString) = @dsp_ccall(
    "readParamFile", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}), 
    dsp.p, param_file)

###############################################################################
# Load problems
###############################################################################

readSmps(dsp::Dsp.Model, filename::AbstractString) = @dsp_ccall(
    "readSmps", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}), 
    dsp.p, convert(Vector{UInt8}, filename))

loadFirstStage(dsp::Dsp.Model, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd) = @dsp_ccall(
    "loadFirstStage", Cvoid,
    (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    dsp.p, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

loadSecondStage(dsp::Dsp.Model, id, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd) = @dsp_ccall(
    "loadSecondStage", Cvoid,
    (Ptr{Cvoid}, Cint, Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    dsp.p, id, probability, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

loadBlockProblem(dsp::Dsp.Model, id, ncols, numels, nrows, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd) = @dsp_ccall(
    "loadBlockProblem", Cvoid, (
    Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt8}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
    dsp.p, id, ncols, nrows, numels, start, index, value, clbd, cubd, ctype, obj, rlbd, rubd)

updateBlocks(dsp::Dsp.Model) = @dsp_ccall("updateBlocks", Cvoid, (Ptr{Cvoid},), dsp.p)


###############################################################################
# solve functions
###############################################################################

for func in [:freeSolver, 
             :solveDe, 
             :solveBd, 
             :solveDd, 
             :solveDw]
    strfunc = string(func)
    @eval begin
        function $func(dsp::Dsp.Model)
            return @dsp_ccall($strfunc, Cvoid, (Ptr{Cvoid},), dsp.p)
        end
    end
end

for func in [:solveBdMpi, :solveDdMpi, :solveDwMpi]
    strfunc = string(func)
    if @isdefined(MPI) #&& MPI.Initialized()
        @eval begin
            function $func(dsp::Dsp.Model, comm)
                return @dsp_ccall($strfunc, Cvoid, (Ptr{Cvoid}, MPI.CComm), dsp.p, MPI.CComm(comm))
	        end
    	end
    else
        @eval begin
            function $func(dsp::Dsp.Model, comm)
                error("MPI package is required to use this function.")
			end
		end
	end
end

###############################################################################
# Get functions
###############################################################################

for (func,rtn) in [(:getNumScenarios, Cint), 
                   (:getTotalNumRows, Cint), 
                   (:getTotalNumCols, Cint), 
                   (:getStatus, Cint), 
                   (:getNumIterations, Cint), 
                   (:getNumNodes, Cint), 
                   (:getWallTime, Cdouble), 
                   (:getPrimalBound, Cdouble), 
                   (:getDualBound, Cdouble),
                   (:getNumCouplingRows, Cint)]
    strfunc = string(func)
    @eval begin
        function $func(dsp::Dsp.Model)
            @assert(dsp.p != C_NULL)
            return @dsp_ccall($strfunc, $rtn, (Ptr{Cvoid},), dsp.p)
        end
    end
end

function getSolution(dsp::Dsp.Model, num::Integer)
    sol = zeros(num)
    if dsp.comm_rank == 0
        @dsp_ccall("getPrimalSolution", Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cdouble}), dsp.p, num, sol)
    end
    return sol
end
getSolution(dsp::Dsp.Model) = getSolution(dsp, getTotalNumCols(dsp))

function getDualSolution(dsp::Dsp.Model, num::Integer)
	sol = zeros(num)
    if dsp.comm_rank == 0
        @dsp_ccall("getDualSolution", Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cdouble}), dsp.p, num, sol)
    end
    return sol
end
getDualSolution(dsp::Dsp.Model) = getDualSolution(dsp, getNumCouplingRows(dsp))

###############################################################################
# Set functions
###############################################################################

setNumberOfScenarios(dsp::Dsp.Model, nscen::Int) = @dsp_ccall(
    "setNumberOfScenarios", Cvoid,
    (Ptr{Cvoid}, Cint), 
    dsp.p, convert(Cint, nscen))

setDimensions(dsp::Dsp.Model, ncols1::Int, nrows1::Int, ncols2::Int, nrows2::Int) = @dsp_ccall(
    "setDimensions", Cvoid,
    (Ptr{Cvoid}, Cint, Cint, Cint, Cint),
    dsp.p, convert(Cint, ncols1), convert(Cint, nrows1), convert(Cint, ncols2), convert(Cint, nrows2))

setIntPtrParam(dsp::Dsp.Model, name::String, n::Int, v::Vector{Int}) = @dsp_ccall(
    "setIntPtrParam", Cvoid, 
    (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}),
    dsp.p, name, convert(Cint, n), convert(Vector{Cint}, v))

end # end of module
