import sys
import math
from input import Data
from gurobipy import *

if len(sys.argv) != 2:
    print("invalid options")
    print("example of use: ./solver.py input.txt")
    exit(1)

# Read input file
data = Data(sys.argv[1])
# nodes = range(data.size)
model = Model("mv-migrate")
# model.setParam(GRB.Param.TimeLimit, 60 * 60 * 10.0) # 10 hour

# Variables
rho = []  # ρ(i) i个vm

for i in range(data.vm_count):
    rho[i] = model.addVar(vtype=GRB.BINARY, name=f"rho_{i}")

# Update model
model.update()

# Add constraints
# TODO

# Objective function
model.setObjective(quicksum(data.vm_bandwidth[i] * data.band_unitcost + data.vm_storage[i] * rho[i] for i in range(data.vm_count)))

# Save LP model and run
model.write('model.rlp')

model.modelSense = GRB.MINIMIZE
model.optimize()

# Print results
if model.status != GRB.OPTIMAL:
    print('No optimum solution found. Status: %i' % (model.status))
else:
    print("Optimal solution found")
    rho_result = [i for i in range(data.vm_count) if rho[i].X > 0.9]
    print("Result = %.4f" % model.objVal)
    print("GAP = %.4f %%" % model.MIPGap)
    print("Time = %.4f seg" % model.Runtime)

    print("ρ(i): ")
    print(rho_result)
