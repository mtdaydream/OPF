from numpy import *
from gurobipy import *
import sys
import time
import copy

casenumber = "toy"
outfile = "CS"+str(casenumber)+"InBdLpFeas"+".csv"
Q = {}
q = {}
c = {}
# Q[0] = array([[1,0],[0,0]])
# Q[1] = array([[0,0],[0,1]])
# q[0] = [-1,0]
# q[1] = [0,1]
# c[0] = 2
# c[1] = 2
Q[0] = array([[1,0],[0,0]])
Q[1] = array([[0,0],[0,1]])
q[0] = [1,-3]
q[1] = [2,-1]
c[0] = -2
c[1] = 4

if casenumber == "toy":
    clist = [h for h in range(2)]

# ceps = 10.0**(-6.0)
# ceps = 0.001
ceps = 1
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

b = [-0.5,3,-0.5,3]
A = array([[-1,0],[1,0],[0,-1],[0,1]])

dimb = len(b)

bqr = []
count = -1
for i in range(dimb):
	for j in range(i,dimb):
		count += 1
		bqr.append([i,j])
		b.append(b[i]*b[j])
		# print(b[i]*b[j])

for k in range(len(clb)):
	b.append(-clb[k])

for k in range(len(clb)):
	b.append(cub[k])	

b.append(b[0])
b.append(-b[0])

dimx = shape(A)[1]
Mrows = len(b)

xqr = []
for i in range(dimx):
	for j in range(i,dimx):
		xqr.append([i,j])

Mcols = shape(A)[1] + len(xqr)
M = zeros((Mrows,Mcols))
for i in range(shape(A)[0]):
	for j in range(shape(A)[1]):
		M[i,j] = A[i,j]

for i in range(len(bqr)):
	bq = bqr[i][0]
	br = bqr[i][1]
	for j in range(shape(A)[1]):
		M[shape(A)[0]+i,j] = b[br]*A[bq,j]+b[bq]*A[br,j]
	for j in range(len(xqr)):
		xi = xqr[j][0]
		xj = xqr[j][1]
		if xi == xj:
			M[shape(A)[0]+i,shape(A)[1]+j] = -A[bq,xi]*A[br,xi]
		else:
			M[shape(A)[0]+i,shape(A)[1]+j] = -A[bq,xj]*A[br,xi]-A[bq,xi]*A[br,xj]

currentrow = shape(A)[0]+i+1

for i in range(len(Q)):
	for j in range(shape(A)[1]):
		M[currentrow+i,j] = -q[i][j]
	for k in range(len(xqr)):
		xi = xqr[k][0]
		xj = xqr[k][1]
		if xi == xj:
			M[currentrow+i,j+k+1] = -Q[i][xi,xj]
		else:
			M[currentrow+i,j+k+1] = -Q[i][xi,xj] - Q[i][xj,xi]


currentrow = currentrow + i+1

for i in range(len(Q)):
	for j in range(shape(A)[1]):
		M[currentrow+i,j] = q[i][j]
	for k in range(len(xqr)):
		xi = xqr[k][0]
		xj = xqr[k][1]
		if xi == xj:
			M[currentrow+i,j+k+1] = Q[i][xi,xj]
		else:
			M[currentrow+i,j+k+1] = Q[i][xi,xj] + Q[i][xj,xi]

currentrow = currentrow + i

for j in range(shape(A)[1]):
	M[currentrow+1,j] = A[0,j]

for j in range(shape(A)[1]):
	M[currentrow+2,j] = -A[0,j]

######## Basic Model Setup
start = time.clock()
## Model parameters and variables
model = Model("InBdFeasModel")
model.setParam('OutputFlag', 0)
### Construct Coefficient matrix M and B
B = b
M = array(M)

### Define model to be feasible by setting objective to 0
# model.setObjective(0.0)

###### Make Model
## Model variables
m = shape(A)[0]
n = shape(A)[1]
lenc = len(clist)
x = dict()
for i in range(shape(M)[1]):
	if i == 0:
		x[i] = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="x"+str(i))
	else:
		x[i] = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, obj=0.0, vtype=GRB.CONTINUOUS, name="x"+str(i))


### Constraints M*X >= B
Constraints = dict()
for k in range(shape(M)[0]-2):
	expr1 = [(M[k,i],x[i]) for i in range(shape(M)[1])]
	Constraints[k] = model.addLConstr(LinExpr(expr1),GRB.LESS_EQUAL, B[k], 'c'+str(k))

expr1 = [(M[shape(M)[0]-2,i],x[i]) for i in range(shape(M)[1])]
Constraints[k+1] = model.addLConstr(LinExpr(expr1),GRB.LESS_EQUAL, B[shape(M)[0]-2], 'c'+str(k+1))

expr1 = [(M[shape(M)[0]-1,i],x[i]) for i in range(shape(M)[1])]
Constraints[k+2] = model.addLConstr(LinExpr(expr1),GRB.LESS_EQUAL, B[shape(M)[0]-1], 'c'+str(k+2))
c_start = shape(A)[0]+len(bqr)
eps = [0.000001,0.00001,0.0001,0.001,0.01,0.1,0.5]

model.update()

while eps:
	success = True
	print("Trying radius, ",c[0]+Constraints[c_start].RHS)
	for i in range(shape(A)[0]):
		for j in range(shape(A)[1]):
			model.chgCoeff(Constraints[shape(M)[0]-2], x[j], A[i,j])
			model.chgCoeff(Constraints[shape(M)[0]-1], x[j], -A[i,j])

		Constraints[shape(M)[0]-2].RHS = b[i]
		Constraints[shape(M)[0]-1].RHS = -b[i]
		model.setObjective(0.0)
		model.optimize()

		try:
			objective = model.objVal
			infeas = False
		except:
			infeas = True

		if infeas:
			print("Row ",i," infeasible. Moving on to row ", i+1 )
		else:
			success = False 
			break
	if success:
		rad = c[0]+Constraints[c_start].RHS
		for i in range(2*len(clb)):
			RHS = Constraints[c_start+i].RHS+ceps
			Constraints[c_start+i].RHS = RHS
		model.update()
	else:
		old_eps = ceps
		if eps:
			ceps = eps.pop()
			for i in range(2*len(clb)):
				origRHS = Constraints[c_start+i].RHS
				Constraints[c_start+i].RHS =  origRHS-old_eps+ceps	
		else:
			for i in range(2*len(clb)):
				origRHS = Constraints[c_start+i].RHS
				Constraints[c_start+i].RHS =  origRHS-old_eps

		model.update()
		

print("InnerBd is ", rad)
