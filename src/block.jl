using Dsp.DspSolverInterface
import JuMP
export @block

type BlockStructure
    parent
    children::Dict{Int,JuMP.Model}
    weight::Dict{Int,Float64}
end

macro block(self, child, id, weight)
    return quote
        if haskey($(esc(self)).ext, :DspBlocks) == false
            $(esc(self)).ext[:DspBlocks] = BlockStructure(
                nothing, Dict{Int,JuMP.Model}(), Dict{Int,Float64}())
        end
        if haskey($(esc(child)).ext, :DspBlocks) == false
            $(esc(child)).ext[:DspBlocks] = BlockStructure(
                $(esc(self)), Dict{Int,JuMP.Model}(), Dict{Int,Float64}())
        else
            $(esc(child)).ext[:DspBlocks].parent = $(esc(self))
        end
        $(esc(self)).ext[:DspBlocks].children[$(esc(id))] = $(esc(child))
        $(esc(self)).ext[:DspBlocks].weight[$(esc(id))]   = $(esc(weight))
        if $(esc(self)).ext[:DspBlocks].parent == nothing
            JuMP.setsolvehook($(esc(self)), DspSolverInterface.dsp_solve)
        end
    end
end