# Kibaek Kim - ANL MCS 2016
# Farmer example from Birge and Louveaux book.

using JuMP, Dsp

NS = 3;                        # number of scenarios
probability = [1/3, 1/3, 1/3]; # probability

CROPS  = 1:3 # set of crops (wheat, corn and sugar beets, resp.)
PURCH  = 1:2 # set of crops to purchase (wheat and corn, resp.)
SELL   = 1:4 # set of crops to sell (wheat, corn, sugar beets under 6K and those over 6K)

Cost     = [150 230 260]   # cost of planting crops
Budget   = 500             # budget capacity
Purchase = [238 210];      # purchase price
Sell     = [170 150 36 10] # selling price
Yield    = [3.0 3.6 24.0;
            2.5 3.0 20.0;
            2.0 2.4 16.0]
Minreq   = [200 240 0]     # minimum crop requirement

# JuMP model
m = Model(NS)

@variable(m, x[i=CROPS,s=1:NS] >= 0, Int)
@objective(m, Min, sum{probability[s] * Cost[i] * x[i,s], i=CROPS, s=1:NS})
@constraint(m, nonant[i=CROPS,s=1], x[i,NS] - x[i,1] == 0)
@constraint(m, nonant[i=CROPS,s=2:NS], x[i,s-1] - x[i,s] == 0)
@constraint(m, const_budget[s=1:NS], sum{x[i,s], i=CROPS} <= Budget)

for s in 1:NS
    blk = Model(m, s, probability[s])

    @variable(blk, y[j=PURCH] >= 0)
    @variable(blk, w[k=SELL] >= 0)
    @objective(blk, Min, sum{Purchase[j] * y[j], j=PURCH} - sum{Sell[k] * w[k], k=SELL})
    @constraint(blk, const_minreq[j=PURCH], Yield[s,j] * x[j,s] + y[j] - w[j] >= Minreq[j])
    @constraint(blk, const_minreq_beets, Yield[s,3] * x[3,s] - w[3] - w[4] >= Minreq[3])
    @constraint(blk, const_aux, w[3] <= 6000)
end
