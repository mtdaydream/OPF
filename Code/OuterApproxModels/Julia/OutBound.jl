using JuMP
using Gurobi

function OutBoundLinX(Q,q,A,b,c,BRad,nmclist, clist)
  # timelim = 300
  x = zeros(length(Q),1)
  dimb = length(b)

  #### Set parameters
  eps = 10.0^(-4)

  n = length(Q)
  Mrows = Int(dimb^2)+dimb
  Mcols = Int(n+n^2)

  ### Construct Verdice Index Dicitonary V and Vb
  V= Dict()
  count = 1
  for i in 1:n
    for j in 1:n
      V[i,j]=count
      count +=1
    end
  end

  Vb= Dict()
  count = 1
  for i in 1:dimb
    for j in 1:dimb
      Vb[i,j]=count
      count +=1
    end
  end

  ### Construct Coefficient matrix M and B
  row = 0
  M=sparse(zeros(Mrows,Mcols))
  B=zeros(1,Mrows)

  for i in 1:dimb
      for j in 1:dimb
          row = Vb[i,j]
          B[1,row] = -b[i]*b[j]
          for k in 1:n
            M[row,k] = -b[j]*A[i,k]-b[i]*A[j,k]
          end
          for l in 1:n
            for m in 1:n
                col = n + V[l,m]
                if m != l
                  M[row,col] = A[j,m]*A[i,l]
                else
                  M[row,col] = A[j,m]*A[i,l]
                end #if
            end #for m
          end #for l
      end #for j
  end #for i

  for i in 1:dimb
      row = Int(dimb^2)+i
      B[1,row] = -b[i]
      for col in 1:n
          M[row,col] = -A[i,col]
      end
  end

  MT = transpose(M)

  #### Construct linear lambda coefficient
  C = Dict()
  for i=1:n
      covec = zeros(1,n)
      row = i
      C[row] = covec
      for k in 1:n
          C[i][k] = q[k][row]
      end
  end

  for L=1:n
    for M=1:n
      row = n+V[L,M]
      covec = zeros(1,n)
      C[row] = covec
      if L!=M
         for k=1:n
             C[row][k] = Q[k][M,L]
         end
      else
        for k=1:n
          C[row][k] = Q[k][M,L]
        end
      end
    end
  end

  ###### Make Model
  # m = Model(solver=GurobiSolver(Presolve=0, TimeLimit = timelim))
  m = Model(solver=GurobiSolver(Presolve=0))
  @variable(m, Y[1:Mrows])
  @variable(m, Lam[1:n])
  for k in nmclist
      @constraint(m, Lam[k] == 0)
  end

  for k in 1:Mcols
    @expression(m, shared1, sum(Y[i]*MT[k,:][i] for i=1:length(MT[k,:])))
    @expression(m, shared2, sum(Lam[i]*C[k][i] for i=1:length(C[k])))
    @constraint(m, shared1 == shared2)
  end

  for k in 1:Mrows
    @constraint(m, Y[k] <= 0)
  end

  #### Only variables in clist should be considered in MIP
  lenc = length(clist)
  @variable(m, LamP[1:lenc])
  @variable(m, LamN[1:lenc])
  @variable(m, category = :Bin, LamS[1:lenc])

  count = 0
  for k in clist
    count += 1
    @constraint(m, LamP[count]-LamN[count] == Lam[k])
    @constraint(m, LamP[count] <= LamS[count])
    @constraint(m, LamP[count] >= 0)
    @constraint(m, LamN[count] <= 1-LamS[count])
    @constraint(m, LamN[count] >= 0)
  end

  @constraint(m, sum(LamP[i]+LamN[i] for i=1:lenc) == 1)

  @objective(m, :Min, sum(B[1,i]*Y[i] for i=1:Mrows))

  status = string(solve(m))
  minobjval = getobjectivevalue(m)

 return minobjval

end

function OutBoundLinXLP(Q,q,A,b,c,BRad,nmclist,clist)
  # timelim = 300
  x = zeros(length(Q),1)
  dimb = length(b)

  #### Set parameters
  eps = 10.0^(-4)

  n = length(Q)
  Mrows = Int(dimb^2)+dimb
  Mcols = Int(n+n^2)

  ### Construct Verdice Index Dicitonary V and Vb
  V= Dict()
  count = 1
  for i in 1:n
    for j in 1:n
      V[i,j]=count
      count +=1
    end
  end

  Vb= Dict()
  count = 1
  for i in 1:dimb
    for j in 1:dimb
      Vb[i,j]=count
      count +=1
    end
  end

  ### Construct Coefficient matrix M and B
  row = 0
  M=sparse(zeros(Mrows,Mcols))
  B=zeros(1,Mrows)

  for i in 1:dimb
      for j in 1:dimb
          row = Vb[i,j]
          B[1,row] = -b[i]*b[j]
          for k in 1:n
            M[row,k] = -b[j]*A[i,k]-b[i]*A[j,k]
          end
          for l in 1:n
            for m in 1:n
                col = n + V[l,m]
                if m != l
                  M[row,col] = A[j,m]*A[i,l]
                else
                  M[row,col] = A[j,m]*A[i,l]
                end #if
            end #for m
          end #for l
      end #for j
  end #for i

  for i in 1:dimb
      row = Int(dimb^2)+i
      B[1,row] = -b[i]
      for col in 1:n
          M[row,col] = -A[i,col]
      end
  end

  MT = transpose(M)

  #### Construct linear lambda coefficient
  C = Dict()
  for i=1:n
      covec = zeros(1,n)
      row = i
      C[row] = covec
      for k in 1:n
          C[i][k] = q[k][row]
      end
  end

  for L=1:n
    for M=1:n
      row = n+V[L,M]
      covec = zeros(1,n)
      C[row] = covec
      if L!=M
         for k=1:n
             C[row][k] = Q[k][M,L]
         end
      else
        for k=1:n
          C[row][k] = Q[k][M,L]
        end
      end
    end
  end

  lambcount= Dict()
  counter = 0
  for i in clist
    counter += 1
    lambcount[i] = counter
  end

  ###### Make Model
  # m = Model(solver=GurobiSolver(Presolve=0, TimeLimit = timelim))
  m = Model(solver=GurobiSolver(Presolve=0))
  @variable(m, Y[1:Mrows])
  @variable(m, Lam[1:length(clist)])

  for k in 1:Mcols
    @expression(m, shared1, sum(Y[i]*MT[k,:][i] for i=1:length(MT[k,:])))
    @expression(m, shared2, sum(Lam[lambcount[i]]*C[k][i] for i in clist))
    @constraint(m, shared1 == shared2)
  end

  for k in 1:Mrows
    @constraint(m, Y[k] <= 0)
  end

  #### Only variables in clist should be considered
  lenc = length(clist)
  @variable(m, LamP[1:lenc])
  @variable(m, LamN[1:lenc])

  for k in clist
    @constraint(m, LamP[lambcount[k]]-LamN[lambcount[k]] == Lam[lambcount[k]])
    @constraint(m, LamP[lambcount[k]] >= 0)
    @constraint(m, LamN[lambcount[k]] >= 0)
  end

  signs= Dict()
  count = 0
  if length(clist) > 5
      for ai in 0:1
          for bi in 0:1
              for ci in 0:1
                  for di in 0:1
                      for ei in 0:1
                          for fi in 0:1
                              for gi in 0:1
                                  for hi in 0:1
                                      for ii in 0:1
                                          for ji in 0:1
                                              sign = [ai,bi,ci,di,ei,fi,gi,hi,ii,ji]
                                              count += 1
                                              signs[count] = sign
                                          end
                                      end
                                  end
                              end
                          end
                      end
                  end
              end
          end
      end

  else
      for ai in 0:1
          for bi in 0:1
              for ci in 0:1
                  for di in 0:1
                      for ei in 0:1
                          sign = [ai,bi,ci,di,ei]
                          count += 1
                          signs[count] = sign
                      end
                  end
              end
          end
      end
  end

  @constraint(m, scon1[k=1:length(clist)], LamP[k] <= signs[1][k])
  @constraint(m, scon2[k=1:length(clist)], LamN[k] <= 1-signs[1][k])

  @constraint(m, sum(LamP[i]+LamN[i] for i=1:lenc) == 1)

  @objective(m, :Min, sum(B[1,i]*Y[i] for i=1:Mrows))

  status = string(solve(m))
  minobjval = getobjectivevalue(m)
  minc = minobjval
  minsign = signs[1]
  for j in 2:length(signs)
    for k in 1:length(clist)
      JuMP.setRHS(scon1[k], signs[j][k])
      JuMP.setRHS(scon2[k], 1-signs[j][k])
    end
    status = string(solve(m))
    minobjval = getobjectivevalue(m)
    if minobjval < minc
      minc = minobjval
      minsign = signs[j]
    end
  end

 return minc

end
