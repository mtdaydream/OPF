from numpy import *
import scipy.io as sio
import sys
import time
import cplex
from time import process_time

'''''''''''''''''''''''''''''''''''
Input/Output Dataset and Initialize
'''''''''''''''''''''''''''''''''''
casenumber = 30
# brad = 0.01
# brad = 0.005
brad = 0.001
file1 = '../../Qscase/Qscase' + str(casenumber) + '.mat'
file2 = '../../Acase/Acase'   + str(casenumber) + '.mat'

out_file = "CS" + str(casenumber) + "PythonMaxZ" + str(brad) + ".csv"

Qs = sio.loadmat(file1, variable_names="Qs")['Qs']
A = sio.loadmat(file2, variable_names="A")['A']
# b = sio.loadmat(file3, variable_names="b")['b']
b = brad * ones(shape(A)[0])

# BigM = max(abs(b))[0]
BigM = 1
dim = int(shape(Qs)[2])
Q = {}
q = {}
c = {}

for k in range(0, dim):
	Q[k] = Qs[1:, 1:, k]
	q[k] = 2*Qs[1:, :1, k]
	c[k] = Qs[0, 0, k]

x = transpose(zeros(shape(Q[0])[1]))

'''''''''''''''''''''
Perform Check
'''''''''''''''''''''

for k in range(0, dim):
	test1 = (transpose(x) * Q[k] * x)[0] + (transpose(x) * q[k])[0] + c[k] - 10.0 ** (-4.0)
	test2 = (transpose(x) * Q[k] * x)[0] + (transpose(x) * q[k])[0] + c[k] + 10.0 ** (-4.0)
	if (0 <= test1).any and (0 >= test2).any:
		print('System not stable')

if casenumber <= 10:
	clist = [x for x in range(5)]
else:
	clist = [x for x in range(10)]

ceps = 10.0 ** (-6.0)
clb = {}
cub = {}

for k in range(dim):
	if k in clist:
		clb[k] = c[k] - ceps
		cub[k] = c[k] + ceps
	else:
		clb[k] = c[k]
		cub[k] = c[k]

nclist = list(range(dim))
for i in clist:
	nclist.remove(i)


# Find maximal c choice
orig_stdout = sys.stdout
f = open(out_file, 'w')
sys.stdout = f
print('Case' + str(casenumber) + 'Linear Constraints using python code, no SDP and uptri vars.')
print('B, CRad, Time, NumVars, NumCons')
sys.stdout = orig_stdout
f.close()


'''''''''''''''''''''
Basic Model Setup
'''''''''''''''''''''

# Model parameters
Gaplim = 0.0
prob = cplex.Cplex()
# prob.parameters.timelimit.set(timelim);
# prob.parameters.mip.limits.solutions.set(Solim);
# prob.parameters.mip.tolerances.integrality.set(0);
prob.parameters.mip.tolerances.mipgap.set(Gaplim)
prob.objective.set_sense(prob.objective.sense.maximize)
# prob.parameters.emphasis.mip.set(1) ## Emphasis set to 1 for Feasibility instead of 0 for balanced


# Model variables
newb = b.copy()
m = shape(A)[0]
n = shape(A)[1]
numvars = int(n * (n + 1) / 2 + n)
numuptri = numvars - n
constraintrowindex = 0
CIndex = {}

# XQ+Xq, first n^2 terms part of XQ, last Xq
xindices = prob.variables.add(names=["x"+str(i) for i in range(numvars + m + 1)])

for x in range(numvars):
	prob.variables.set_lower_bounds(x, -cplex.infinity)

for x in range(m):
	prob.variables.set_types(numvars+x, prob.variables.type.binary)

prob.variables.set_lower_bounds(numvars+m, -cplex.infinity)

# Create Upper Triangular dictionary
UpTriDic = {}
counter = 0
for i in range(n):
	for j in range(i, n):
		UpTriDic[i, j] = counter
		UpTriDic[j, i] = counter
		counter += 1


# Constraint Ax<=b
my_rhs = newb
prob.linear_constraints.add(rhs=my_rhs)
for L in range(m):
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["axb", L] = x
	prob.linear_constraints.set_senses(x, "L")
	for R in range(n):
		y = numuptri + R
		v = A[L, R]
		prob.linear_constraints.set_coefficients(x, y, v)


# Constraint -Axb^T-b(Ax)^T+AXA^T >= -bb
for L in range(m):
	for K in range(m):
		x = constraintrowindex
		constraintrowindex += 1
		CIndex["axa", L, K] = x
		my_rhs = [-1 * newb[L] * newb[K]]
		prob.linear_constraints.add(rhs=my_rhs)
		prob.linear_constraints.set_senses(x, "G")
		for r in range(n):
			y = numuptri + r
			v = -(A[L, r] * newb[K] + newb[L] * A[K, r])
			prob.linear_constraints.set_coefficients(x, y, v)
			for s in range(r, n):
				y = UpTriDic[r, s]
				if r != s:
					v = A[K, s] * A[L, r] + A[K, r] * A[L, s]
				else:
					v = A[K, s] * A[L, r]
				prob.linear_constraints.set_coefficients(x, y, v)

# Tr(QX) + q^tx <=/>= -clb/-cub
for k in clist:
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["qxcl", k] = x
	my_rhs = [-1 * clb[k]]
	prob.linear_constraints.add(rhs=my_rhs)
	prob.linear_constraints.set_senses(x, "L")
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["qxcu", k] = x
	my_rhs = [-1 * cub[k]]
	prob.linear_constraints.add(rhs=my_rhs)
	prob.linear_constraints.set_senses(x, "G")
	for r in range(n):
		y = numuptri + r
		v = q[k][r, 0]
		prob.linear_constraints.set_coefficients(x-1, y, v)
		prob.linear_constraints.set_coefficients(x, y, v)
		for s in range(r, n):
			y = UpTriDic[r, s]
			if r != s:
				v = Q[k][r, s] + Q[k][s, r]
			else:
				v = Q[k][s, s]
			prob.linear_constraints.set_coefficients(x-1, y, v)
			prob.linear_constraints.set_coefficients(x, y, v)

	
for k in nclist:
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["qxcl", k] = x
	my_rhs = [-1 * clb[k]]
	prob.linear_constraints.add(rhs=my_rhs)
	prob.linear_constraints.set_senses(x, "L")
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["qxcu", k] = x
	my_rhs = [-1 * cub[k]]
	prob.linear_constraints.add(rhs=my_rhs)
	prob.linear_constraints.set_senses(x, "G")
	for r in range(n):
		y = numuptri + r
		v = q[k][r, 0]
		prob.linear_constraints.set_coefficients(x-1, y, v)
		prob.linear_constraints.set_coefficients(x, y, v)
		for s in range(r, n):
			y = UpTriDic[r, s]
			if r != s:
				v = Q[k][r, s] + Q[k][s, r]
			else:
				v = Q[k][s, s]
			prob.linear_constraints.set_coefficients(x-1, y, v)
			prob.linear_constraints.set_coefficients(x, y, v)

# Constraint z <= A[i,:]x-b_i+BigM*(1-d_i) (d_i binary)
for i in range(m):
	x = constraintrowindex
	constraintrowindex += 1
	CIndex["BigM", i] = x
	my_rhs = [b[i] - BigM]
	prob.linear_constraints.add(rhs=my_rhs)
	prob.linear_constraints.set_senses(x, "G")
	prob.linear_constraints.set_coefficients(x, numvars+m, -1)
	prob.linear_constraints.set_coefficients(x, numvars+i, -BigM)
	for j in range(n):
		y = numuptri + j
		v = A[i, j]
		prob.linear_constraints.set_coefficients(x, y, v)

# Constraint sum d_i =1
x = constraintrowindex
constraintrowindex += 1
CIndex["SumBin"] = x
my_rhs = [1]
prob.linear_constraints.add(rhs=my_rhs)
prob.linear_constraints.set_senses(x, "E")
for i in range(m):
	y = numvars + i
	v = 1
	prob.linear_constraints.set_coefficients(x, y, v)

prob.objective.set_linear(numvars+m, 1)

# Loop Parameters
eps = 10.0 ** (-4.0)
eps2 = eps / 10.0
eps3 = 0.01   # Incremental change in C bounds

if brad <= eps:
	eps = brad / 10.0
	eps2 = eps / 10.0
	
maxcub = cub.copy()
minclb = clb.copy()
testcub = cub.copy()
testclb = clb.copy()

prob.solve()
z = prob.solution.get_objective_value()

if z >= -0.00001:
	print('Initial ceps of %d is too large' % ceps)
	quit()

for k in clist:
	testclb[k] -= eps3
	testcub[k] += eps3

# Reset constraint Tr(QX) + q^tx <=/>= -clb/-cub
for k in clist:
	x = CIndex["qxcl", k]
	my_rhs = -1 * testclb[k]
	prob.linear_constraints.set_rhs(x, my_rhs)

	x = CIndex["qxcu", k]
	my_rhs = -1 * testcub[k]
	prob.linear_constraints.set_rhs(x, my_rhs)

NonConvergentRuns = 0
start = process_time()

while NonConvergentRuns <= 2:
	if testcub[0] > 1:
		print('max cub reached')
		quit()

	prob.solve()
	z = prob.solution.get_objective_value()
	# print z, -0.00001
	# print z+0.00001
	# time.sleep(1)

	# if str(type(z)) != "<type 'float'>":
	if isinstance(z,float) is False:
		print('Model stopped \n******************************************************\n')
		print('Found NaN \n******************************************************')
		quit()

	elif z < -0.00001:
		maxcub = testcub.copy()
		minclb = testclb.copy()

		for k in clist:
			testclb[k] -= eps3
			testcub[k] += eps3

		# Reset constraint Tr(QX) + q^tx <=/>= -clb/-cub
		for k in clist:
			x = CIndex["qxcl", k]
			my_rhs = -1 * testclb[k]
			prob.linear_constraints.set_rhs(x, my_rhs)

			x = CIndex["qxcu", k]
			my_rhs = -1 * testcub[k]
			prob.linear_constraints.set_rhs(x, my_rhs)

	else:
		NonConvergentRuns += 1
		for k in clist:
			testclb[k] += eps3
			testcub[k] -= eps3

		eps3 = eps3/10.0
		for k in clist:
			testclb[k] -= eps3
			testcub[k] += eps3

		# Reset constraint Tr(QX) + q^tx <=/>= -clb/-cub
		for k in clist:
			x = CIndex["qxcl", k]
			my_rhs = -1*testclb[k]
			prob.linear_constraints.set_rhs(x, my_rhs)

			x = CIndex["qxcu", k]
			my_rhs = -1*testcub[k]
			prob.linear_constraints.set_rhs(x, my_rhs)

# End find Max C Rad, print to file
etime = process_time() - start
icrad = maxcub[clist[0]]
numv = prob.variables.get_num()
numc = prob.linear_constraints.get_num()
print('max c rad is ', icrad)
orig_stdout = sys.stdout
f = open(out_file, 'a')
sys.stdout = f
print(brad, ',', icrad, ',', etime, ",", numv, ",", numc)
sys.stdout = orig_stdout
f.close()