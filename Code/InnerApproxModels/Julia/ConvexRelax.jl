using JuMP
# using Gurobi
using Mosek

# Might nor work with Julia 0.6 -> Need to fix vectorized constraints

function CvxRlxLinX(Q,q,clb,cub,A,b,Cord)
  # m = Model(solver=GurobiSolver(Presolve=0))
  m = Model(solver=MosekSolver())
  sq =length(Q)
  sA = size(A)[1]
  B = b*(b')
  @variable(m, Xq[1:sq])
  @variable(m, XQ[1:sq,1:sq], Symmetric)
  @SDconstraint(m, XQ >= zeros(sq,sq))

  ### Ax <= b
  @constraint(m, A*Xq .<= b)

  ### bb-Axb^T-b(Ax)^T+AXA^T >= 0
  exp3 = @expression(m, A*XQ*(A'))
  B = b*b'
  exp2 = @expression(m, A*Xq*b')
  exp1 = @expression(m, b*(A*Xq)')
  z = zeros(sA)
  @constraint(m, B-exp1-exp2+exp3 .>= z)

  ## Tr(QX)+q'x+clb <=0   & Tr(QX)+q'x+cub >=0
  for k in 1:sq
      @constraint(m, trace(Q[k]*XQ) + q[k]'*Xq+clb[k] .<= 0)
      @constraint(m, trace(Q[k]*XQ) + q[k]'*Xq+cub[k] .>= 0)
  end

  ### Max Ax at dim Cord
  objexp = slice(A, Cord, :)
  @objective(m, :Max, (objexp'*Xq)[1])
  status = string(solve(m))

  return getobjectivevalue(m), status
end
