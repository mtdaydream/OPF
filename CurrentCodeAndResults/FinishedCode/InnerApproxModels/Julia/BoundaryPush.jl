include("ConvexRelax.jl") # Implements the model in JuMP - builds the relaxation

function CBdPushLinX(Q,q,c,A,b,clb,cub,clist)
  eps = 10.0^(-4)   ###epsilon bound for convergence between bounds
  eps2 = 10.0^(-5)  ###epsilon bound for insufficient bound change
  if b[1] <= max(eps,eps2)
      eps = b[1]/10.0
      eps2 = eps/10.0
  end
  eps2 = 0.0 ### For absolute case
  eps3 = 0.001   ###epsilon bound for incremental change in c's
  dim = length(b)
  NonConvergentruns = 0
  origcub = copy(cub)
  origclb = copy(clb)
  minclb = copy(clb)
  maxcub = copy(cub)
  ob = copy(b)

  while NonConvergentruns < 2
      for k in clist
        cub[k] += eps3
        clb[k] -= eps3
      end

      origcub = copy(cub)
      origclb = copy(clb)

      TestStatus = 1
      DimCheck = ones(dim,1)
      b = copy(ob)

      while TestStatus > 0
        oldb = copy(b)
        @sync @parallel for d in 1:dim

          if DimCheck[d,1] == 1
             Nb, status = CvxRlxLinX(Q,q,clb,cub,A,b,d)

             if length(Nb) > 1
              TestStatus = -1
              println("Status is ", status)
              println("******************************************************")
              println("Found NaN")
              println("******************************************************")
              sleep(5)
              quit()

            elseif b[d]-Nb > 0
              if Nb <= eps2
                b[d] = eps2
                DimCheck[d,1] = 0
              else
                b[d] = Nb
              end  ## End eps if check
            end ## End Nb if check
          end ## End dimension if check
        end ## End dimension for loop

        if maximum(DimCheck) == 0 || minimum(abs(ob-b)) > eps
          TestStatus = 0
          maxcub = origcub
          minclb = origclb

       elseif maximum(oldb-b) <= eps2
          TestStatus = 0
          NonConvergentruns += 1
          if NonConvergentruns <= 2
            for k in clist
                cub[k] -= eps3
                clb[k] += eps3
            end
            eps3 = eps3/10.0
          end
       end  ##### End TestStatus Check
   end  #### End while TestStatus While loop
 end   ##### End while loop NonConvergentruns
 return minclb, maxcub
end
