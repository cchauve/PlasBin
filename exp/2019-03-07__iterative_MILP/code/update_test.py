from __future__ import division
from gurobipy import *


w = [2,3,4,5]
k_list = [1,3]
k = [0,1,0,1]

m = Model("Update")

x = {}
for i in range(4):
	x[i] = m.addVar(vtype=GRB.BINARY, name='var_'+str(i))

expr = LinExpr()
for i in x:
	expr.addTerms(0.5*w[i],x[i])

m.setObjective(expr, GRB.MAXIMIZE)

constr_1 = LinExpr()
constr_2 = LinExpr()
for i in x:
	constr_1.addTerms(1,x[i])
	constr_2.addTerms(k[i],x[i])

m.addConstr(constr_1 <= 2, "constr_"+str(1))
m.addConstr(constr_2 >= 1, "constr_"+str(2))


#m.optimize()

while len(k_list) > 0:
	m.optimize()
	print(k_list)
	constr_new = LinExpr()
	for i in x:
		if x[i].x > 0:
			print(x[i].varName, x[i].x, 0.5*w[i])

			w[i] = w[i]/2

			if w[i] <= 2:
				if k[i] == 1:
					k[i] = 0
					k_list.remove(i)
		constr_new.addTerms(k[i],x[i])				
					
	m.addConstr(constr_new >= 1, "constr_new")

	expr = LinExpr()
	for i in x:
		expr.addTerms(0.5*w[i],x[i])

	m.setObjective(expr, GRB.MAXIMIZE)

