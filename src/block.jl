import JuMP

type BlockStructure
    parent
    children::Dict{Int,JuMP.Model}
    weight::Dict{Int,Float64}
end