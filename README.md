# DSP.jl
[![Build Status](https://travis-ci.org/kibaekkim/DSP.jl.svg?branch=master)](https://travis-ci.org/kibaekkim/DSP.jl)
[![codecov](https://codecov.io/gh/kibaekkim/Dsp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kibaekkim/Dsp.jl)

DSP.jl is an interface to a parallel decomposition mixed-integer programming solver [DSP](https://github.com/Argonne-National-Laboratory/DSP). 
This package allows users to define block structures in optimization model written in [StructJuMP](https://github.com/StructJuMP/StructJuMP.jl) 
and solve the block-structured problem using the parallle solver ``DSP``.

## Intallation

> **_NOTE:_** You need to install solver [DSP](https://github.com/Argonne-National-Laboratory/DSP) first. This package provides an interface only.

Due to the dependency on `StructJuMP.jl`, we recommend to instantiate the package with `Manifest.toml`.
To do so, please `clone` and `instantiate` the repo.

```
git clone git@github.com:kibaekkim/DSP.jl.git
cd DSP.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Examples

This is one simple example of `stochastic` form.
Please find more examples in `./examples` particularly for `block` form.

```julia
using MPI
using StructJuMP
using DSP

# Comment out this line if you want to run in serial
MPI.Init()

# Initialize DSP.jl with the communicator.
DSP.parallelize(MPI.COMM_WORLD)

xi = [[7,7] [11,11] [13,13]]

# create StructuredModel with number of scenarios
m = StructuredModel(num_scenarios = 3)

@variable(m, 0 <= x[i=1:2] <= 5, Int)
@objective(m, Min, -1.5 * x[1] - 4 * x[2])
for s = 1:3
    # create a StructuredModel linked to m with id s and probability 1/3
    blk = StructuredModel(parent = m, id = s, prob = 1/3)
    @variable(blk, y[j=1:4], Bin)
    @objective(blk, Min, -16 * y[1] + 19 * y[2] + 23 * y[3] + 28 * y[4])
    @constraint(blk, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= xi[1,s] - x[1])
    @constraint(blk, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= xi[2,s] - x[2])
end

status = optimize!(m, 
    is_stochastic = true, # Needs to indicate that the model is a stochastic program.
    solve_type = DSP.Dual, # see instances(DSP.Methods) for other methods
)

# Comment out this line if you want to run in serial
MPI.Finalize()
```

## Acknowledgements
This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.
