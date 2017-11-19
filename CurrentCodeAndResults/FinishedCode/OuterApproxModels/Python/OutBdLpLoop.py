from numpy import *
import scipy.io as sio
import sys
import time
import cplex
import copy
import scipy as sci 

casenumber = 22
# brad = 0.001
# brad = 0.005
brad = 0.01
file1 = "C:/Users/Ben Rapone/Documents/Research/OPF/PNNL/DJRepo/ccsi-robustoptfeas/Ben/Code/src/Qscase"+str(casenumber)+".mat";
file2 = "C:/Users/Ben Rapone/Documents/Research/OPF/PNNL/DJRepo/ccsi-robustoptfeas/Ben/Code/src/Acase"+str(casenumber)+".mat";
file3 = "C:/Users/Ben Rapone/Documents/Research/OPF/PNNL/DJRepo/ccsi-robustoptfeas/Ben/Code/src/bcase"+str(casenumber)+".mat";

# file1 = "Qscase"+str(casenumber)+".mat";
# file2 = "Acase"+str(casenumber)+".mat";
# file3 = "bcase"+str(casenumber)+".mat";


outfile = "CS"+str(casenumber)+"OutBdLpLoopBR"+str(brad)+".csv"

Qs = sio.loadmat(file1,variable_names="Qs")['Qs']
A = sio.loadmat(file2,variable_names="A")['A']
# b = sio.loadmat(file3,variable_names="b")['b']
b= brad*ones(shape(A)[0])


dim = int(shape(Qs)[2])
Q = {}
q = {}
c = {}

for k in range(0,dim):
    Q[k] = Qs[1:,1:,k]
    q[k] = 2*Qs[1:,:1,k]
    c[k] = Qs[0,0,k]

x = transpose(matrix(zeros((shape(Q[0])[1]))))

# #### Perform Check

for k in range(0,dim):
  if 0 <= (transpose(x)*Q[k]*x)[0]+(transpose(x)*q[k])[0]+c[k]-10.0**(-4.0) and 0 >= (transpose(x)*Q[k]*x)[0]+(transpose(x)*q[k])[0]+c[k]+10.0**(-4.0):
    print "System not stable"

if casenumber <= 100:
    clist = [x for x in range(5)]
else:
    clist = [x for x in range(10)]

ceps = 10.0**(-6.0)
clb = {}
cub = {}

for k in range(dim):
  if k in clist:
      clb[k] = c[k]-ceps
      cub[k] = c[k]+ceps
  else:
      clb[k] = c[k]
      cub[k] = c[k]

nclist = range(dim)
for i in clist:
	nclist.remove(i)

##### Find maximal c choice

orig_stdout = sys.stdout
f = file(outfile,'w')
sys.stdout = f
print "Case "+ str(casenumber)+" Linear Constraints using python code with no SDP and uptri vars. Outer Bound estimate."
print "Brad",",", "Minc",",", "Time"
sys.stdout = orig_stdout
f.close()

######## Basic Model Setup
start = time.clock()
## Model parameters
Gaplim = 0.0
prob = cplex.Cplex();
# prob.parameters.timelimit.set(timelim);
# prob.parameters.mip.limits.solutions.set(Solim);
# prob.parameters.mip.tolerances.integrality.set(0);
# prob.parameters.mip.tolerances.mipgap.set(Gaplim);
prob.objective.set_sense(prob.objective.sense.minimize);

#### Set parameters
# timelim = 300
x = zeros((len(Q),1))
dimb = len(b)
eps = 10.0**(-4.0)

n = len(Q)
Mrows = int(dimb**2)+dimb
Mcols = int(n+n**2)

### Construct Verdice Index Dicitonary V and Vb
V= dict()
count = 0
for i in range(0,n): 
	for j in range(0,n):
		V[i,j]=count
		count +=1

Vb= dict()
count = 0
for i in range(0,dimb):
	for j in range(0,dimb):
		Vb[i,j]=count
		count +=1

### Construct Coefficient matrix M and B
row = 0
M=sci.sparse.lil_matrix((Mrows,Mcols))
B=zeros((1,Mrows))

for i in range(0,dimb):
	for j in range(0,dimb):
		row = Vb[i,j]
		B[0,row] = -b[i]*b[j]
		for k in range(n):
			M[row,k] = -b[j]*A[i,k]-b[i]*A[j,k]

		for l in range(n):
			for m in range(n):
				col = n + V[l,m]
				if m != l:
					M[row,col] = A[j,m]*A[i,l]
				else:
					M[row,col] = A[j,m]*A[i,l]

for i in range(dimb):
	row = int(dimb**2)+i
	B[0,row] = -b[i]
	for col in range(n):
		M[row,col] = -A[i,col]

MT = transpose(M)

#### Construct linear lambda coefficient
C = dict()
for row in range(n):
	covec = zeros(n)
	C[row] = covec
	for k in range(n):
		C[row][k] = q[k][row]

for L in range(n):
	for M in range(n):
		row = n+V[L,M]
		covec = zeros(n)
		C[row] = covec
		if L!=M:
			for k in range(n):
				C[row][k] = Q[k][M,L]
		else:
			for k in range(n):
				C[row][k] = Q[k][M,L]

signs= dict()
count = 0
if len(clist) > 5:
	for ai in range(2):
		for bi in range(2):
			for ci in range(2):
				for di in range(2):
					for ei in range(2):
						for fi in range(2):
							for gi in range(2):
								for hi in range(2):
									for ii in range(2):
										for ji in range(2):
											sign = [ai,bi,ci,di,ei,fi,gi,hi,ii,ji]
											signs[count] = sign
											count += 1
else:
	for ai in range(2):
		for bi in range(2):
			for ci in range(2):
				for di in range(2):
					for ei in range(2):
						sign = [ai,bi,ci,di,ei]
						signs[count] = sign
						count += 1


###### Make Model
## Model variables
m = shape(A)[0]
n = shape(A)[1]
lenc = len(clist)
numvars = Mrows+3*lenc #### Y variables are first Mrow vars, lam vars are next lenc vars, Positive lam vars are next lenc vars, Negative lam vars are remaining lenc vars
constraintrowindex = 0
CIndex = {}
yindices = prob.variables.add(names = ["Y"+str(i) for i in range(Mrows)])  
for y in range(Mrows):
    prob.variables.set_lower_bounds(y, -cplex.infinity)
    prob.variables.set_upper_bounds(y, 0)

lamindices = prob.variables.add(names = ["Lam"+str(i) for i in range(lenc)])
for y in range(lenc):
    prob.variables.set_lower_bounds(Mrows+y, -cplex.infinity)

Plamindices = prob.variables.add(names = ["PLam"+str(i) for i in range(lenc)])  
Nlamindices = prob.variables.add(names = ["NLam"+str(i) for i in range(lenc)])  

### Constraints MT[k,:]*Y = Lam*C 
my_rhs = zeros(Mcols)
prob.linear_constraints.add(rhs=my_rhs)
for k in range(Mcols):
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["MTY-LC",k] = x
	prob.linear_constraints.set_senses(x,"E")
	nz = sci.sparse.lil_matrix.nonzero(MT[k,:])[1]
	for y in nz:
		v = MT[k,y]
		prob.linear_constraints.set_coefficients(x,int(y),v)
	for i in clist:
		y = Mrows+i
		v = -C[k][i]
		prob.linear_constraints.set_coefficients(x,y,v)

#### Lambda Constraints Lam = LamP-LamN
my_rhs = zeros(lenc)
prob.linear_constraints.add(rhs=my_rhs)
for k in clist:
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["Lams",k] = x
	y = Mrows+k
	v = -1
	prob.linear_constraints.set_coefficients(x,y,v)
	y = Mrows+lenc+k
	v = 1
	prob.linear_constraints.set_coefficients(x,y,v)
	y = Mrows+2*lenc+k
	v = -1
	prob.linear_constraints.set_coefficients(x,y,v)
	prob.linear_constraints.set_senses(x,"E")


#### Sum LamP+LamN ==1
my_rhs = ones(1)
prob.linear_constraints.add(rhs=my_rhs)
x = constraintrowindex
constraintrowindex += 1
CIndex["SumLamPN",k] = x
prob.linear_constraints.set_senses(x,"E")
for k in clist:
	y = Mrows+lenc+k
	v = 1
	prob.linear_constraints.set_coefficients(x,y,v)
	y = Mrows+2*lenc+k
	v = 1
	prob.linear_constraints.set_coefficients(x,y,v)

#### Sign Constraints
for k in clist:
	my_rhs = [signs[0][k]]
	prob.linear_constraints.add(rhs=my_rhs)
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["SignsP",k] = x
	y = Mrows+lenc+k
	v = 1
	prob.linear_constraints.set_coefficients(x,y,v)
	prob.linear_constraints.set_senses(x,"L")

for k in clist:
	my_rhs = [1-signs[0][k]]
	prob.linear_constraints.add(rhs=my_rhs)
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["SignsL",k] = x
	y = Mrows+2*lenc+k
	v = 1
	prob.linear_constraints.set_coefficients(x,y,v)
	prob.linear_constraints.set_senses(x,"L")

for i in range(Mrows):
	prob.objective.set_linear(i,B[0,i])

prob.solve()
print prob.solution.get_objective_value()
# print prob.linear_constraints.get_rows(78)

minobjval = prob.solution.get_objective_value()
minc = minobjval
minsign = signs[0]

for j in range(1,len(signs)):
	for k in clist:
		x = CIndex["SignsP", k]
		my_rhs = signs[j][k]
		prob.linear_constraints.set_rhs(x,my_rhs)
		x = CIndex["SignsL", k]
		my_rhs = 1-signs[j][k]
		prob.linear_constraints.set_rhs(x,my_rhs)
	prob.solve()
	minobjval = prob.solution.get_objective_value()
	# print minobjval
	if minobjval < minc:
		minc = minobjval
etime = time.clock()-start
# print "Minc is ", minc

orig_stdout = sys.stdout
f = file(outfile,'a')
sys.stdout = f
print brad,",", minc, ",", etime
sys.stdout = orig_stdout
f.close()






