using DSP
using Test

const dsp = DSP.dspenv

@testset "Initializing DSP" begin
    @test dsp.p != C_NULL
    @test length(dsp.numRows) == 0
    @test length(dsp.numCols) == 0
    @test isnan(dsp.primVal)
    @test isnan(dsp.dualVal)
    @test length(dsp.colVal) == 0
    @test length(dsp.rowVal) == 0
    @test dsp.nblocks == 0
    @test dsp.block_ids == []
    @test dsp.is_stochastic == false
    @test dsp.solve_type == DSP.Dual
    @test isnothing(dsp.comm)
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
end

@testset "Setting options" begin
    @testset "param:" begin
        DSP.setoptions!(Dict(:param => "params.txt"))
    end
    @testset "is_stochastic: $i" for i in [true, false]
        DSP.setoptions!(Dict(:is_stochastic => i))
        @test dsp.is_stochastic == i
    end
    @testset "solve_type: $t" for t in instances(DSP.Methods)
        DSP.setoptions!(Dict(:solve_type => t))
        @test dsp.solve_type == t
    end
end

@testset "Farmer example: stochastic form" begin
    include("farmer_stoc.jl")
    @testset "Parent model" begin
        @test length(m.variables) == 3
        @test length(m.constraints) == 1
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSP.get_model_data(m)
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
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSP.get_model_data(subm)
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
        @testset "load_problem!" begin
            DSP.load_problem!(m)
            @test DSP.getNumSubproblems(dsp) == 3
            # @test DSP.getTotalNumRows(dsp) == 13
            @test DSP.getTotalNumCols(dsp) == 21
        end
        @testset "Methods: $j" for j in instances(DSP.Methods)
            @testset "solve!" begin
                dsp.solve_type = j
                DSP.solve!()
                @test dsp.status == 3000
                @test DSP.termination_status(m) == MOI.OPTIMAL
            end
            @testset "post_solve!" begin
                dsp.solve_type = j
                DSP.post_solve!()

                if dsp.solve_type in [DSP.Dual, DSP.ExtensiveForm]
                    @test isapprox(objective_value(m), -108390.)
                else
                    @test objective_value(m) >= -108390.
                end
                @test isapprox(dual_objective_value(m), -108390.)

                primsol = value()
                dualsol = dual()
                if dsp.solve_type == DSP.Legacy
                    for k = 0:3
                        @test primsol[k] != []
                    end
                    @test dualsol != []
                else
                    @test isapprox(primsol[0], [170.0, 80.0, 250.0])
                    if dsp.solve_type != DSP.Benders
                        @test isapprox(primsol[1], [0.0, 0.0, 310.0, 48.0, 6000.0, 0.0])
                        @test isapprox(primsol[2], [0.0, 0.0, 225.0, 0.0, 5000.0, 0.0])
                        @test isapprox(primsol[3], [0.0, 48.0, 140.0, 0.0, 4000.0, 0.0])
                    end
                    @test dualsol == []
                    @test isapprox(value(x[1]), 170.0)
                    @test isapprox(value(x[2]), 80.0)
                    @test isapprox(value(x[3]), 250.0)
                end
            end
            @testset "freeSolver" begin
                DSP.freeSolver(dsp)
            end
        end
        @testset "freeModel" begin
            DSP.freeModel(dsp)
            @test dsp.p != C_NULL
            @test length(dsp.numRows) == 0
            @test length(dsp.numCols) == 0
            @test isnan(dsp.primVal)
            @test isnan(dsp.dualVal)
            @test length(dsp.colVal) == 0
            @test length(dsp.rowVal) == 0
            @test dsp.nblocks == 0
            @test dsp.block_ids == []
            @test dsp.is_stochastic == false
            @test dsp.solve_type == DSP.Dual
        end
    end
end

@testset "Farmer example: block form" begin
    include("farmer_block.jl")
    @testset "Parent model" begin
        @test length(m.variables) == 27
        @test length(m.constraints) == 6
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSP.get_model_data(m)
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
        start, index, value, rlbd, rubd, obj, clbd, cubd, ctype, cname = DSP.get_model_data(subm)
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
        @testset "load_problem!" begin
            DSP.load_problem!(m)
            @test DSP.getNumSubproblems(dsp) == 3
            # @test DSP.getTotalNumRows(dsp) == 18
            @test DSP.getTotalNumCols(dsp) == 27
        end
        @testset "Methods: $j" for j in [DSP.ExtensiveForm, DSP.Dual]
            @testset "solve!" begin
                dsp.solve_type = j
                DSP.solve!()
                @test dsp.status == 3000
            end
            @testset "post_solve!" begin
                dsp.solve_type = j
                DSP.post_solve!()

                @test isapprox(objective_value(m), -108390.)
                @test isapprox(dual_objective_value(m), -108390.)

                primsol = value()
                dualsol = dual()
                @test isapprox(primsol[0], [
                    170.0, 170.0, 170.0, 
                    80.0, 80.0, 80.0, 
                    250.0, 250.0, 250.0, 
                    0.0, 0.0, 0.0, 
                    0.0, 0.0, 48.0, 
                    310.0, 225.0, 140.0, 
                    48.0, 0.0, 0.0, 
                    6000.0, 5000.0, 4000.0, 
                    0.0, 0.0, 0.0])
                for s = 1:3
                    @test primsol[s] == []
                end
                @test dualsol == []
                @test isapprox(value(x[1,1]), 170.0)
                @test isapprox(value(x[1,2]), 170.0)
                @test isapprox(value(x[1,3]), 170.0)
                @test isapprox(value(x[2,1]), 80.0)
                @test isapprox(value(x[2,2]), 80.0)
                @test isapprox(value(x[2,3]), 80.0)
                @test isapprox(value(x[3,1]), 250.0)
                @test isapprox(value(x[3,2]), 250.0)
                @test isapprox(value(x[3,3]), 250.0)
                @test isapprox(value(y[1,1]), 0.0)
                @test isapprox(value(y[1,2]), 0.0)
                @test isapprox(value(y[1,3]), 0.0)
                @test isapprox(value(y[2,1]), 0.0)
                @test isapprox(value(y[2,2]), 0.0)
                @test isapprox(value(y[2,3]), 48.0)
                @test isapprox(value(w[1,1]), 310.0)
                @test isapprox(value(w[1,2]), 225.0)
                @test isapprox(value(w[1,3]), 140.0)
                @test isapprox(value(w[2,1]), 48.0)
                @test isapprox(value(w[2,2]), 0.0)
                @test isapprox(value(w[2,3]), 0.0)
                @test isapprox(value(w[3,1]), 6000.0)
                @test isapprox(value(w[3,2]), 5000.0)
                @test isapprox(value(w[3,3]), 4000.0)
                @test isapprox(value(w[4,1]), 0.0)
                @test isapprox(value(w[4,2]), 0.0)
                @test isapprox(value(w[4,3]), 0.0)
            end
            @testset "freeSolver" begin
                DSP.freeSolver(dsp)
            end
        end
        @testset "freeModel" begin
            DSP.freeModel(dsp)
            @test dsp.p != C_NULL
            @test length(dsp.numRows) == 0
            @test length(dsp.numCols) == 0
            @test isnan(dsp.primVal)
            @test isnan(dsp.dualVal)
            @test length(dsp.colVal) == 0
            @test length(dsp.rowVal) == 0
            @test dsp.nblocks == 0
            @test dsp.block_ids == []
            @test dsp.is_stochastic == false
            @test dsp.solve_type == DSP.Dual
        end
    end
end

@testset "Freeing DSP" begin
    DSP.freeEnv(dsp)
    @test dsp.p == C_NULL
    @test length(dsp.numRows) == 0
    @test length(dsp.numCols) == 0
    @test isnan(dsp.primVal)
    @test isnan(dsp.dualVal)
    @test length(dsp.colVal) == 0
    @test length(dsp.rowVal) == 0
    @test dsp.nblocks == 0
    @test dsp.block_ids == []
    @test dsp.is_stochastic == false
    @test dsp.solve_type == DSP.Dual
    @test isnothing(dsp.comm)
    @test dsp.comm_size == 1
    @test dsp.comm_rank == 0
end
