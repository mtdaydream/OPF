#This is the script that you run to run the algorithm and get inner and outer bounds on c

include("BoundaryPush.jl") # Pushes the c bounds to limits
include("MakeSystem.jl") # Reads or generates input data
include("OutBound.jl") # Finds outer c bounds limits

casenumber = 57
file = string("C:/Users/Ben Rapone/Documents/Research/PNNL/DJRepo/ccsi-robustoptfeas/Ben/Code/src/Qscase",casenumber,".mat");
file2 = string("C:/Users/Ben Rapone/Documents/Research/PNNL/DJRepo/ccsi-robustoptfeas/Ben/Code/src/Acase",casenumber,".mat");
file3 = string("C:/Users/Ben Rapone/Documents/Research/PNNL/DJRepo/ccsi-robustoptfeas/Ben/Code/src/bcase",casenumber,".mat");

  # x, Q, q, c, ub, lb, cub, clb = GenSys(dim)
x, Q, q, c, cub, clb, A, b, clist = ImportSys(file,file2,file3)
dim = length(Q)


nmclist = [x for x=1:length(c)]
nmclist = setdiff(nmclist,clist)
ocub = copy(cub)
oclb = copy(clb)

b = b - 0.1*ones(length(b))
# bsteps = [0.001,0.002,0.003, 0.004,0.005,0.006,0.007,0.008,0.009,0.01]
bsteps = [0.001]

#### Inner Bound
# GraphInfo = Dict()
# Outfile = string("CS",casenumber,"LinConIn.csv")
# open(Outfile, "w") do f
# write(f, string("Case",casenumber," with linear constraints SDP"), "\n")
# close(f)
# end
# for i in 1:Int(length(bsteps))
#     brad = bsteps[i]
#     newb = copy(b) + brad*ones(length(b))
#     tic()
#     minclb, maxcub = CBdPushLinX(Q,q,c,A,newb,clb,cub,clist)
#     etime = toc()
#     icrad = maxcub[clist[1]][1]-c[clist[1]][1]
#     GraphInfo[i]=[brad,icrad,etime]
#     open(Outfile, "a") do f
#     write(f, "$(GraphInfo[i][1])")
#     for j in 2:3
#         write(f, ",$(GraphInfo[i][j])")
#     end
#     write(f, "\n")
#     close(f)
#     end
# end

# #### Outer Bound
# GraphInfo2 = Dict()
# Outfile = string("CS",casenumber,"LinConOut.csv")
# open(Outfile, "w") do f
# write(f, string("Case",casenumber," with linear constraints no SDP"), "\n")
# close(f)
# end
# for i in 1:Int(length(bsteps))
#     brad = bsteps[i]
#     newb = copy(b) + brad*ones(length(b))
#     tic()
#     maxc = OutBoundLinX(Q,q,A,newb,c,brad,nmclist,clist)
#     etime = toc()
#     GraphInfo2[i]=[brad,maxc,etime]
#     open(Outfile, "a") do f
#     write(f, "$(GraphInfo2[i][1])")
#     for j in 2:3
#         write(f, ",$(GraphInfo2[i][j])")
#     end
#     write(f, "\n")
#     close(f)
#     end
# end

# #### Outer Bound LP Loop
GraphInfo2 = Dict()
Outfile = string("CS",casenumber,"LinConOutLP.csv")
open(Outfile, "w") do f
write(f, string("Case",casenumber," with linear constraints no SDP"), "\n")
close(f)
end
for i in 1:Int(length(bsteps))
    brad = bsteps[i]
    newb = copy(b) + brad*ones(length(b))
    tic()
    maxc = OutBoundLinXLP(Q,q,A,newb,c,brad,nmclist,clist)
    etime = toc()
    GraphInfo2[i]=[brad,maxc,etime]
    open(Outfile, "a") do f
    write(f, "$(GraphInfo2[i][1])")
    for j in 2:3
        write(f, ",$(GraphInfo2[i][j])")
    end
    write(f, "\n")
    close(f)
    end
end
