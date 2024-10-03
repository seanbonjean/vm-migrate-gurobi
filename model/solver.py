import sys
import math
from input import Data
from gurobipy import *

if len(sys.argv) != 2:
    print("invalid options")
    print("example of use: ./solver.py input.txt")
    exit(1)

'''
def find_lesshop(vm_num):
    cluster_num = d.vm_mycluster[vm_num]
    same_cluster = []
    # 寻找同集群的vm
    for i in range(d.vm_count):
        if d.vm_mycluster[i] == cluster_num and i != vm_num:
            same_cluster.append(i)
    # 如果同集群为空（无同集群vm），则返回0（说明无需vm间通信）
    if same_cluster == []:
        return 0
    # 寻找同集群vm中跳数最少的vm
    lesshop = math.inf
    for i in same_cluster:
        hop = d.hop_matrix[(d.vm_mypm[vm_num], d.vm_mypm[i])]
        if hop < lesshop:
            lesshop = hop
            target_vm = i
            target_pm = d.vm_mypm[i]
    return (lesshop, target_vm, target_pm)
'''


def cal_hops(vm_num):
    cluster_num = d.vm_mycluster[vm_num]
    same_cluster = []
    # 寻找同集群的vm（自身除外）
    for i in range(d.vm_count):
        if d.vm_mycluster[i] == cluster_num and i != vm_num:
            same_cluster.append(i)
    # 如果同集群为空（无同集群vm），则返回0（说明无需vm间通信）
    if same_cluster == []:
        return 0
    # 对同集群除自身外所有vm，计算跳数并累加
    total_hop = 0
    for i in same_cluster:
        hop = d.hop_matrix[(omiga[vm_num], d.vm_mypm[i])]
        total_hop += hop
    return total_hop / 2


def cal_Rmi(Umi):
    return 1+math.log(1-(1-1/math.e)*Umi)


# Read input file
d = Data(sys.argv[1])
model = Model("mv-migrate")
# model.setParam(GRB.Param.TimeLimit, 60 * 60 * 10.0) # 10 hour

# Variables
rho = {}  # ρ(i) i个vm是否迁移
omiga = {}  # ω(i) vm迁移到目标pm的位置

for i in range(d.vm_count):
    rho[i] = model.addVar(vtype=GRB.BINARY, name=f"rho_{i}")

for i in range(d.vm_count):
    omiga[i] = model.addVar(
        vtype=GRB.INTEGER, name=f"omiga_{i}", lb=1, ub=d.pm_count)

# Update model
model.update()

# Add constraints
# 大M法
M = 10000000
b = [model.addVar(vtype=GRB.BINARY, name=f"b_{i}") for i in range(d.vm_count)]
for i in range(d.vm_count):
    model.addConstr((rho[i] == 0) >> (omiga[i] == d.vm_mypm[i]))
    # model.addConstr((rho[i] == 1) >> (omiga[i] != d.vm_mypm[i]))
    model.addConstr((rho[i] == 1) >> (
        omiga[i] <= d.vm_mypm[i] - 1 + M * (1 - b[i])))
    # 当 b[i] = 1 时，强制 omega[i] <= d.vm_mypm[i] - 1
    model.addConstr((rho[i] == 1) >> (omiga[i] >= d.vm_mypm[i] + 1 - M * b[i]))
    # 当 b[i] = 0 时，强制 omega[i] >= d.vm_mypm[i] + 1

for i in range(1, d.pm_count + 1):
    pm_core_sum = quicksum(d.vm_CPUcore[j] for j in range(d.vm_count) if omiga[j] == i)
    Uc = quicksum(d.vm_CPUuti[j] * d.vm_CPUcore[j] for j in range(d.vm_count) if omiga[j] == i) / pm_core_sum
    Uo = pm_core_sum / d.pm_CPUcore[i-1]
    Ur = quicksum(d.vm_storage[j] for j in range(d.vm_count) if omiga[j] == i) / d.pm_storage[i-1]
    Ub = quicksum(d.vm_bw[j] for j in range(d.vm_count) if omiga[j] == i) / d.pm_bw[i-1]
    model.addConstr(cal_Rmi(Uc) >= d.elasticity_lb)
    model.addConstr(cal_Rmi(Uo) >= d.elasticity_lb)
    model.addConstr(cal_Rmi(Ur) >= d.elasticity_lb)
    model.addConstr(cal_Rmi(Ub) >= d.elasticity_lb)
    # model.addConstr(pm_core_sum <= d.pm_CPUcore[i-1])
    model.addConstr(Uc <= 1)
    model.addConstr(Uo <= 1)
    model.addConstr(Ur <= 1)
    model.addConstr(Ub <= 1)
    model.addConstr(Uc >= 0)
    model.addConstr(Uo >= 0)
    model.addConstr(Ur >= 0)
    model.addConstr(Ub >= 0)

# Objective function
model.setObjective(quicksum(d.cluster_innerpakage_size[d.vm_mycluster[i]] * cal_hops(
    i) * d.band_unitcost + d.vm_storage[i] * rho[i] for i in range(d.vm_count)))

# Save LP model and run
model.write('model.rlp')

model.modelSense = GRB.MINIMIZE
model.optimize()

# Print results
if model.status != GRB.OPTIMAL:
    print('No optimum solution found. Status: %i' % (model.status))
else:
    print("Optimal solution found")
    rho_result = [1 if rho[i].X > 0.9 else 0 for i in range(d.vm_count)]
    omiga_result = [omiga[i].X for i in range(d.vm_count)]
    print("Result = %.4f" % model.objVal)
    print("GAP = %.4f %%" % model.MIPGap)
    print("Time = %.4f seg" % model.Runtime)

    print("ρ(i): ")
    print(rho_result)
    print("ω(i): ")
    print(omiga_result)
