# Dsp.jl

Dsp.jl is an interface to a parallel decomposition mixed-integer programming solver [DSP](https://github.com/Argonne-National-Laboratory/DSP). This package allows users to define block structures in optimization model written in [JuMP](https://github.com/JuliaOpt/JuMP.jl) and solve the block-structured problem using the parallle solver ``DSP``.

# Intallation

```julia
Pkg.clone("https://github.com/kibaekkim/Dsp.jl")
```

# Example

DSP can read and solve model from JuMP:

```julia
using JuMP, Dsp, MPI

# Comment out this line if you want to run in serial
MPI.Init()

xi = [[7,7] [11,11] [13,13]]

# create JuMP.Model with number of blocks
m = Model(3)

@variable(m, 0 <= x[i=1:2] <= 5, Int)
@objective(m, Min, -1.5 * x[1] - 4 * x[2])
for s = 1:3
    # create a JuMP.Model block linked to m with id s and probability 1/3
    blk = Model(m, s, 1/3)
    @variable(blk, y[j=1:4], Bin)
    @objective(blk, Min, -16 * y[1] + 19 * y[2] + 23 * y[3] + 28 * y[4])
    @constraint(blk, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= xi[1,s] - x[1])
    @constraint(blk, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= xi[2,s] - x[2])
end

solve_types = [:Dual, :Benders, :Extensive]
solve(m, solve_type = solve_types[1], param = "myparam.txt")

getobjectivevalue(m)

# Comment out this line if you want to run in serial
MPI.Finalize()

```

or, it can also read and solve model from SMPS files:

```julia
using Dsp

# Assumming we have mysmps.cor, mysmps.tim, and mysmps.sto
readSmps("mysmps")

solve_types = [:Dual, :Benders, :Extensive]
optimize(solve_type = solve_types[1], param = "myparam.txt")

println(getprimobjval())
println(getdualobjval())
```
