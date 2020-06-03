using Dsp
using Test

const dsp = Dsp.dsp_model

@testset "Initializing Model" begin
    @test dsp.p != C_NULL
    @test dsp.solve_type == :Dual
    @test dsp.numRows == 0
    @test dsp.numCols == 0
    @test isnan(dsp.primVal)
    @test isnan(dsp.dualVal)
    @test dsp.colVal == []
    @test dsp.rowVal == []
    @test dsp.comm == nothing
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
end

@testset "Setting options" begin
    for t in [:Dual, :Benders]
        Dsp.setoptions(dsp, Dict(:solve_type => t))
        @test dsp.solve_type == t
    end
end

# include("../examples/farmer.jl")
# Dsp.optimize!(m)
