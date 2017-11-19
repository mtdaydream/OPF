using MAT

# Reads mat file into the format used in JuMP code
function ImportSys(fileAdd,fileAdd2,fileAdd3)
  file = matopen(fileAdd)
  Qs = read(file, "Qs")
  close(file)

  file = matopen(fileAdd2)
  A = read(file, "A")
  close(file)

  file = matopen(fileAdd3)
  b = read(file, "b")
  b= 0.1*ones(length(b))
  close(file)

  len = Int(size(Qs[:,:,1])[1])
  Q = Dict()
  q = Dict()
  c = Dict()

  for k in 1:len-1
    Q[k] = slice(Qs[:,:,k], 2:len, 2:len)
    q[k] = 2*slice(Qs[:,:,k], 2:len, 1)
    c[k] = Qs[1,1,k]
  end
  x = zeros(len-1,1)


 #### Perform Check
 # println("Ax-b is ", A*x-b)
  for k in 1:len-1
    if 0 <= (x'*Q[k]*x)[1]+(x'q[k])[1]+c[k]-10.0^(-4) && 0 >= (x'*Q[k]*x)[1]+(x'q[k])[1]+c[k]-10.0^(-4)
      println("System not stable")
      quit()
    end
  end

  dim = len-1

  if len <= 17
      clist = [x for x=1:5]
  else
      clist = [x for x=1:10]
  end

  eps = 10.0^(-4)
  # lb = x-eps*ones(dim,1)
  # ub = x+eps*ones(dim,1)
  ceps = 10.0^(-6)
  clb = Dict()
  cub = Dict()

    for k in 1:dim
    clb[k] = c[k]
    cub[k] = c[k]
    end

    for k in clist
    clb[k] = c[k]-ceps
    cub[k] = c[k]+ceps
    end

  return x, Q, q, c, cub, clb, A, b, clist
end

# Generates a random system with appropriate dimensions
function GenSys(dim)
  ceps = 10.0^(-6)
  eps = 10.0^(-4)
  x = randn(dim,1)
  lb = x-eps*ones(dim,1)
  ub = x+eps*ones(dim,1)

  Qi = randn(dim,dim)
  qi = randn(dim,1)
  ci = -(0.5*transpose(x)*Qi*x+transpose(qi)*x)


  c = Dict([(1,ci)])
  Q = Dict([(1,Qi)])
  q = Dict([(1,qi)])
  clb = Dict([(1,ci-ceps)])
  cub = Dict([(1,ci+ceps)])
  # clb = Dict()
  # cub = Dict()
  # clb[1] = ci[1]-ceps
  # cub[1] = ci[1]+ceps

  for k in 2:dim
    Qi = randn(dim,dim)
    qi = randn(dim,1)
    ci = -(0.5*transpose(x)*Qi*x+transpose(qi)*x)
    c[k] = ci
    clb[k] = ci-ceps
    cub[k] = ci+ceps
    q[k] = qi
    Q[k] = Qi
  end

  return x, Q, q, c, ub, lb, cub, clb
end
