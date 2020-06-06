using Dsp
using Test

const dsp = Dsp.dsp_model

@testset "Initializing DSP" begin
    @test dsp.p != C_NULL
    @test dsp.solve_type == Dsp.Dual
    @test dsp.numRows == 0
    @test dsp.numCols == 0
    @test isnan(dsp.primVal)
    @test isnan(dsp.dualVal)
    @test dsp.colVal == []
    @test dsp.rowVal == []
    @test dsp.comm == nothing
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
    @test dsp.block_ids == []
    @test dsp.is_stochastic == false
end

@testset "Setting options" begin
    @testset "param:" begin
        Dsp.setoptions(dsp, Dict(:param => "params.txt"))
    end
    @testset "is_stochastic: $i" for i in [true, false]
        Dsp.setoptions(dsp, Dict(:is_stochastic => i))
        @test dsp.is_stochastic == i
    end
    @testset "solve_type: $t" for t in instances(Dsp.Methods)
        Dsp.setoptions(dsp, Dict(:solve_type => t))
        @test dsp.solve_type == t
    end
end
#=
@testset "Farmer example: stochastic form" begin
    include("../examples/farmer_stoc.jl")
    @testset "Parent model" begin
        @test length(m.variables) == 3
        @test length(m.constraints) == 1
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = Dsp.get_model_data(m)
        @test length(start) == length(m.constraints) + 1
        @test start[end] == length(index)
        @test start == [0, 3]
        @test index == [0, 1, 2]
        @test value == [1., 1., 1.]
        @test rlbd == [-Inf]
        @test rubd == [Budget]
        @test obj == Cost
        @test clbd == zeros(3)
        @test cubd == zeros(3) .+ Inf
        @test ctype == "III"
    end
    @testset "Child model $i" for (i,subm) in m.children
        @test i > 0
        @test m.probability[i] == probability[i]
        @test length(subm.variables) == 6
        @test length(subm.constraints) == 4
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = Dsp.get_model_data(subm)
        @test length(start) == length(subm.constraints) + 1
        @test start[end] == length(index)
        @test start == [0, 3, 6, 9, 10]
        @test index == [0, 3, 5, 1, 4, 6, 2, 7, 8, 7]
        @test value == [Yield[i,1], 1., -1., Yield[i,2], 1., -1., Yield[i,3], -1., -1., 1.]
        @test rlbd == [Minreq; -Inf]
        @test rubd == [Inf, Inf, Inf, 6000]
        @test obj == [Purchase; -Sell]
        @test clbd == zeros(6)
        @test cubd == zeros(6) .+ Inf
        @test ctype == "CCCCCC"
    end
    @testset "optimize!" begin
        dsp.is_stochastic = true
        @testset "loadProblem" begin
            Dsp.loadProblem(dsp, m)
        end
        @testset "Methods: $j" for j in instances(Dsp.Methods)
            @testset "solve" begin
                dsp.solve_type = j
                Dsp.solve(dsp)
            end
            @show Dsp.DspCInterface.getStatus(dsp)
            @show Dsp.DspCInterface.getPrimalBound(dsp)
            @show Dsp.DspCInterface.getDualBound(dsp)
            @testset "freeSolver" begin
                Dsp.freeSolver(dsp)
            end
        end
        @testset "freeModel" begin
            Dsp.freeModel(dsp)
            @test dsp.p != C_NULL
            @test dsp.solve_type == Dsp.Dual
            @test dsp.numRows == 0
            @test dsp.numCols == 0
            @test isnan(dsp.primVal)
            @test isnan(dsp.dualVal)
            @test dsp.colVal == []
            @test dsp.rowVal == []
            @test dsp.block_ids == []
            @test dsp.is_stochastic == false
        end
    end
end
=#
@testset "Farmer example: block form" begin
    include("../examples/farmer_block.jl")
    @testset "Parent model" begin
        @test length(m.variables) == 27
        @test length(m.constraints) == 6
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = Dsp.get_model_data(m)
        @test length(start) == length(m.constraints) + 1
        @test start[end] == length(index)
        @test start == [0, 2, 4, 6, 8, 10, 12]
        @test index == [0, 1, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8]
        @test value == [1., -1., 1., -1., 1., -1., 1., -1., 1., -1., 1., -1.]
        @test rlbd == zeros(6)
        @test rubd == zeros(6)
        @test obj == [
            probability[1]*Cost[1], probability[2]*Cost[1], probability[3]*Cost[1],
            probability[1]*Cost[2], probability[2]*Cost[2], probability[3]*Cost[2],
            probability[1]*Cost[3], probability[2]*Cost[3], probability[3]*Cost[3],
            probability[1]*Purchase[1], probability[2]*Purchase[1], probability[3]*Purchase[1],
            probability[1]*Purchase[2], probability[2]*Purchase[2], probability[3]*Purchase[2],
            -probability[1]*Sell[1], -probability[2]*Sell[1], -probability[3]*Sell[1],
            -probability[1]*Sell[2], -probability[2]*Sell[2], -probability[3]*Sell[2],
            -probability[1]*Sell[3], -probability[2]*Sell[3], -probability[3]*Sell[3],
            -probability[1]*Sell[4], -probability[2]*Sell[4], -probability[3]*Sell[4]]
        @test clbd == zeros(27)
        @test cubd == zeros(27) .+ Inf
        @test ctype == "IIIIIIIIICCCCCCCCCCCCCCCCCC"
    end
    @testset "Child model $i" for (i,subm) in m.children
        @test i > 0
        @test length(subm.variables) == 0
        @test length(subm.constraints) == 5
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = Dsp.get_model_data(subm)
        @test length(start) == length(subm.constraints) + 1
        @test start[end] == length(index)
        @test start == [0, 3, 6, 9, 12, 13]
        @test index == [0, 3, 6,  0, 9, 15,  3, 12, 18,  6, 21, 24,  21] .+ (i-1)
        @test value == [1., 1., 1.,  Yield[i,1], 1., -1.,  Yield[i,2], 1., -1.,  Yield[i,3], -1., -1., 1.]
        @test rlbd == [-Inf; Minreq; -Inf]
        @test rubd == [Budget, Inf, Inf, Inf, 6000]
        @test obj == []
        @test clbd == []
        @test cubd == []
        @test ctype == ""
    end
    @testset "optimize!" begin
        dsp.is_stochastic = false
        @testset "loadProblem" begin
            Dsp.loadProblem(dsp, m)
        end
        @testset "Methods: $j" for j in [Dsp.ExtensiveForm, Dsp.Dual]
            @testset "solve" begin
                dsp.solve_type = j
                Dsp.solve(dsp)
            end
            @show Dsp.DspCInterface.getStatus(dsp)
            @show Dsp.DspCInterface.getPrimalBound(dsp)
            @show Dsp.DspCInterface.getDualBound(dsp)
            @testset "freeSolver" begin
                Dsp.freeSolver(dsp)
            end
        end
        @testset "solve" begin
            dsp.solve_type = Dsp.ExtensiveForm
            Dsp.solve(dsp)
        end
        @testset "freeSolver" begin
            Dsp.freeSolver(dsp)
        end
        @testset "freeModel" begin
            Dsp.freeModel(dsp)
            @test dsp.p != C_NULL
            @test dsp.solve_type == Dsp.Dual
            @test dsp.numRows == 0
            @test dsp.numCols == 0
            @test isnan(dsp.primVal)
            @test isnan(dsp.dualVal)
            @test dsp.colVal == []
            @test dsp.rowVal == []
            @test dsp.block_ids == []
            @test dsp.is_stochastic == false
        end
    end
end

@testset "Freeing DSP" begin
    Dsp.freeEnv(dsp)
    @test dsp.p == C_NULL
    @test dsp.solve_type == Dsp.Dual
    @test dsp.numRows == 0
    @test dsp.numCols == 0
    @test isnan(dsp.primVal)
    @test isnan(dsp.dualVal)
    @test dsp.colVal == []
    @test dsp.rowVal == []
    @test dsp.comm == nothing
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
    @test dsp.block_ids == []
    @test dsp.is_stochastic == false
end
