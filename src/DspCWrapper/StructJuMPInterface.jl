module IpoptInterfaceSerial

using StructJuMP, JuMP, StructJuMPSolverInterface
using Ipopt

import MathProgBase

type StructJuMPModel <: ModelInterface
    model::JuMP.Model 

    write_solution::Function
    get_x0::Function
    numvars::Function
    numcons::Function
    get_bounds::Function

    function StructJuMPModel(model)
        instance = new(model)

        instance.write_solution = function(x)
            @assert length(x) == getTotalNumVars(m)
            m = instance.model
            idx = 1
            for i = 0:num_scenarios(m)
                mm = getModel(m,i)
                for j = 1:getNumVars(m,i)
                    setvalue(Variable(mm,j),x[idx])
                    idx += 1
                end
            end
        end

        instance.get_x0 = function(x)
            @assert length(x) == getTotalNumVars(m)
            m = instance.model
            idx = 1
            for i = 0:num_scenarios(m)
                mm = getModel(m,i)
                for j = 1:getNumVars(m,i)
                    v_j = getvalue(Variable(mm,j))
                    x[idx] = isnan(v_j)? 1.0:v_j
                    idx += 1
                end
            end
            # @show x
            return x
        end
        instance.numvars = function()
            return getTotalNumVars(instance.model)
        end

        instance.numcons = function()
            return getTotalNumCons(instance.model)
        end

        return instance  
    end
end


function structJuMPSolve(model; suppress_warmings=false, kwargs...)
    # @show typeof(model)
    nm = StructJuMPModel(model)
    n = getTotalNumVars(model)
    m = getTotalNumCons(model)
    nele_jac = nm.nele_jac()
    nele_hess = nm.nele_hess()

    # @show x_L, x_U
    # @show g_L, g_U
    # @show n,m
    # @show nele_jac,nele_hess

    prob = createProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                         nm.eval_f, nm.eval_g, nm.eval_grad_f, nm.eval_jac_g, nm.eval_h)
    # setProblemScaling(prob,1.0)
    nm.get_x0(prob.x)
    status = solveProblem(prob)
    nm.write_solution(prob.x)
    return status
end

KnownSolvers["Ipopt"] = IpoptInterfaceSerial.structJuMPSolve

end  #end module
