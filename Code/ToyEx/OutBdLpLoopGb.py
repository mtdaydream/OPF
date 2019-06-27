from numpy import *
from gurobipy import *
import sys
import time
import copy

casenumber = "toy"

# outfile = "CS"+str(casenumber)+"OutBdLpLoopBR"+str(brad)+".csv"
Q = {}
q = {}
c = {}
Q[0] = array([[1,0],[0,0]])
Q[1] = array([[0,0],[0,1]])
q[0] = [-1,0]
q[1] = [0,1]
c[0] = 2
c[1] = 2

b = [-0.5,3,-0.5,3]
newA = array([[-1,0],[1,0],[0,-1],[0,1]])

b.append(1)
b.append(-1)
Arows = len(b)
Acols = shape(newA)[1]+1
dimb = len(b)-2

bqr = []
count = -1
for i in range(dimb):
	for j in range(i,dimb):
		count += 1
		bqr.append([i,j])
		b.append(b[i]*b[j])

A = zeros((Arows,Acols))
for i in range(shape(newA)[0]):
	for j in range(shape(newA)[1]):
		A[i,j]=newA[i,j]

A[-2,-1] = 1
A[-1,-1] = -1

dimx = shape(A)[1]
Mrows = len(b)

xqr = []
for i in range(dimx-1):
	for j in range(i,dimx-1):
		xqr.append([i,j])


Mcols = shape(A)[1] + len(xqr)
M = zeros((Mrows,Mcols))
for i in range(shape(A)[0]):
	for j in range(shape(A)[1]):
		M[i,j] = A[i,j]

for i in range(len(bqr)):
	bq = bqr[i][0]
	br = bqr[i][1]
	for j in range(shape(A)[1]-1):
		M[shape(A)[0]+i,j] = b[br]*A[bq,j]+b[bq]*A[br,j]
	for j in range(len(xqr)):
		xi = xqr[j][0]
		xj = xqr[j][1]
		if xi == xj:
			M[shape(A)[0]+i,shape(A)[1]+j] = -A[bq,xi]*A[br,xi]
		else:
			M[shape(A)[0]+i,shape(A)[1]+j] = -A[bq,xj]*A[br,xi]-A[bq,xi]*A[br,xj]

M = array(M)
if casenumber == "toy":
    clist = [h for h in range(2)]

ceps = 10.0**(-6.0)
clb = {}
cub = {}
dim = 2
for k in range(dim):
  if k in clist:
      clb[k] = c[k]-ceps
      cub[k] = c[k]+ceps
  else:
      clb[k] = c[k]
      cub[k] = c[k]

######## Basic Model Setup
start = time.clock()
## Model parameters and variables
model = Model("OutBdModel")

### Construct Coefficient matrix M and B
B = b
MT = transpose(M)
	
###### Make Model
## Model variables
m = shape(A)[0]
n = shape(A)[1]-1
lenc = len(clist)
x = dict()

for i in range(Mrows):
	x[i] = model.addVar(lb=0.0, ub=GRB.INFINITY, obj=b[i], vtype=GRB.CONTINUOUS, name="x"+str(i))

Lam = model.addVars(lenc, lb=-GRB.INFINITY, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS)
PLam = model.addVars(lenc,lb=0, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS)

### Constraints MT*Y - LamC = 0 
for k in range(shape(MT)[0]):
	expr1 = [(MT[k,i],x[i]) for i in range(shape(MT)[1])]
	if k < n:
		expr2 = [(-q[i][k],Lam[i]) for i in range(lenc)]
		model.addLConstr(LinExpr(expr1+expr2),GRB.GREATER_EQUAL, 0.0, 'c'+str(k))
	elif k==n:
		expr2 = [(c[i],Lam[i]) for i in range(lenc)]
		model.addLConstr(LinExpr(expr1+expr2),GRB.GREATER_EQUAL, 0.0, 'c'+str(k))
	else:
		xi = xqr[k-n-1][0]
		xj = xqr[k-n-1][1]
		if xi == xj:
			expr2 = [(-Q[i][xi,xj],Lam[i]) for i in range(lenc)]
			model.addLConstr(LinExpr(expr1+expr2),GRB.GREATER_EQUAL, 0.0, 'c'+str(k))
		else:
			expr2 = [(-Q[i][xi,xj]-Q[i][xj,xi],Lam[i]) for i in range(lenc)]
			model.addLConstr(LinExpr(expr1+expr2),GRB.GREATER_EQUAL, 0.0, 'c'+str(k))

#### PLam == abs(Lam)
for i in range(lenc):
	model.addConstr(PLam[i] == abs_(Lam[i]), name = 'c'+str(k+i))

#### Sum PLam 
expr = [(1,PLam[i]) for i in range(lenc)]
model.addLConstr(LinExpr(expr),GRB.EQUAL, 1.0, 'c'+str(k+i+1))

model.write("model_form.rlp")

model.optimize()

model.write("model_solution_Out.sol")
for v in model.getVars():
	print('%s %g' % (v.varName, v.x))

print('Obj: %g' % model.objVal)
