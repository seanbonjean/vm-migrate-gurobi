import sys
import math
from input import Data
from gurobipy import *
from func import calE_star

'''
if len(sys.argv) != 2:
    print("invalid options")
    print("example of use: ./solver.py input.txt")
    exit(1)
'''

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
        # hop = d.hop_matrix[(omiga[vm_num], d.vm_mypm[i])]
        hop = quicksum(gamma[vm_num, j] * d.hop_matrix[(j-1, d.vm_mypm[i]-1)]
                       for j in range(1, d.pm_count + 1))
        total_hop += hop
    return total_hop / 2  # 两vm间共享同一信道，而此处会被计算2次，因此除以2


def cal_Rmi(Umi):
    return 1+math.log(1-(1-1/math.e)*Umi)


def cal_Rmi_reverse(x):
    return (1 - math.exp(x - 1)) / (1 - 1 / math.e)


# Read input file
# d = Data(sys.argv[1])
d = Data("data/test.txt")
model = Model("mv-migrate")
# model.setParam(GRB.Param.TimeLimit, 60 * 60 * 10.0) # 10 hour

# Variables
rho = {}  # ρ(i) i个vm是否迁移
omiga = {}  # ω(i) vm迁移到目标pm的位置

for i in range(d.vm_count):
    rho[i] = model.addVar(vtype=GRB.BINARY, name=f"rho_{i}")
'''
for i in range(d.vm_count):
    omiga[i] = model.addVar(
        vtype=GRB.INTEGER, name=f"omiga_{i}", lb=1, ub=d.pm_count)
'''

# Update model
model.update()

# Add constraints
# !!!下面涉及了两个变换，用(1)表示变换1的注释，用(2)表示变换2的注释，用(3)表示变换3的注释
# 变换(1)用于解决gurobipy.GurobiError: Constraint has no bool value (are you trying "lb <= expr <= ub"?)
# 变换(2)用于解决gurobipy.GurobiError: Divisor must be a constant
# 变换(3)用于解决TypeError: must be real number, not gurobipy.LinExpr

# (1) 定义二进制变量 gamma[j, i] 表示 VM j 是否分配到 PM i
gamma = {}
for j in range(d.vm_count):
    for i in range(1, d.pm_count + 1):
        gamma[j, i] = model.addVar(vtype=GRB.BINARY, name=f"gamma_{j}_{i}")
model.update()

# (2) 引入新的变量 Uc_var 来表示 Uc
Uc_var = {}
for i in range(1, d.pm_count + 1):
    Uc_var[i] = model.addVar(vtype=GRB.CONTINUOUS, name=f"Uc_var_{i}")
model.update()

# (1) 添加约束：只有当 omiga[j] == i 时，gamma[j, i] == 1
for j in range(d.vm_count):
    model.addConstr(quicksum(gamma[j, i] for i in range(
        1, d.pm_count + 1)) == 1, name=f"omiga_assign_{j}")

Uc = {}
Uo = {}
Ur = {}
Ub = {}
for i in range(1, d.pm_count + 1):
    # (1) pm_core_sum = quicksum(d.vm_CPUcore[j] for j in range(d.vm_count) if omiga[j] == i)
    # (1) Uc = quicksum(d.vm_CPUuti[j] * d.vm_CPUcore[j] for j in range(d.vm_count) if omiga[j] == i) / pm_core_sum
    # (1) Uo = pm_core_sum / d.pm_CPUcore[i-1]
    # (1) Ur = quicksum(d.vm_storage[j] for j in range(d.vm_count) if omiga[j] == i) / d.pm_storage[i-1]
    # (1) Ub = quicksum(d.vm_bw[j] for j in range(d.vm_count) if omiga[j] == i) / d.pm_bw[i-1]
    # (1) 报错的约束中，用 gamma[j, i] 替代 omiga[j] == i
    # (1) 使用 gamma[j, i] 来选择 VM 分配给 PM i 的项
    pm_core_sum = quicksum(d.vm_CPUcore[j] * gamma[j, i]
                           for j in range(d.vm_count))
    # (2) Uc = quicksum(d.vm_CPUuti[j] * d.vm_CPUcore[j] * gamma[j, i] for j in range(d.vm_count)) / pm_core_sum
    model.addConstr(quicksum(d.vm_CPUuti[j] * d.vm_CPUcore[j] * gamma[j, i]
                    for j in range(d.vm_count)) == Uc_var[i] * pm_core_sum)
    Uo[i] = pm_core_sum / d.pm_CPUcore[i-1]
    Ur[i] = quicksum(d.vm_storage[j] * gamma[j, i]
                     for j in range(d.vm_count)) / d.pm_storage[i-1]
    Ub[i] = quicksum(d.vm_bw[j] * gamma[j, i]
                     for j in range(d.vm_count)) / d.pm_bw[i-1]

    # (2) model.addConstr(cal_Rmi(Uc) >= d.elasticity_lb)
    # (3) model.addConstr(cal_Rmi(Uc_var[i]) >= d.elasticity_lb)
    # (3) model.addConstr(cal_Rmi(Uo) >= d.elasticity_lb)
    # (3) model.addConstr(cal_Rmi(Ur) >= d.elasticity_lb)
    # (3) model.addConstr(cal_Rmi(Ub) >= d.elasticity_lb)
    R = cal_Rmi_reverse(d.elasticity_lb)
    model.addConstr(Uc_var[i] <= R)
    model.addConstr(Uo[i] <= R)
    model.addConstr(Ur[i] <= R)
    model.addConstr(Ub[i] <= R)
    # model.addConstr(pm_core_sum <= d.pm_CPUcore[i-1])
    # (2) model.addConstr(Uc <= 1)
    model.addConstr(Uc_var[i] <= 1)
    model.addConstr(Uo[i] <= 1)
    model.addConstr(Ur[i] <= 1)
    model.addConstr(Ub[i] <= 1)
    # (2) model.addConstr(Uc >= 0)
    model.addConstr(Uc_var[i] >= 0)
    model.addConstr(Uo[i] >= 0)
    model.addConstr(Ur[i] >= 0)
    model.addConstr(Ub[i] >= 0)

'''
# 大M法，用于解决gurobipy.GurobiError: Inequality constraints not supported报错
# !!!如果pm的数量级比较大，需要把M的值设置得更大
M = 10000000
b = [model.addVar(vtype=GRB.BINARY, name=f"b_{i}") for i in range(d.vm_count)]
model.update()
for i in range(d.vm_count):
    model.addConstr((rho[i] == 0) >> (omiga[i] == d.vm_mypm[i]))
    # model.addConstr((rho[i] == 1) >> (omiga[i] != d.vm_mypm[i]))
    model.addConstr((rho[i] == 1) >> (
        omiga[i] <= d.vm_mypm[i] - 1 + M * (1 - b[i])))
    # 当 b[i] = 1 时，强制 omega[i] <= d.vm_mypm[i] - 1
    model.addConstr((rho[i] == 1) >> (omiga[i] >= d.vm_mypm[i] + 1 - M * b[i]))
    # 当 b[i] = 0 时，强制 omega[i] >= d.vm_mypm[i] + 1
'''
for i in range(d.vm_count):
    model.addConstr((rho[i] == 0) >> (gamma[i, d.vm_mypm[i]] == 1))
    model.addConstr((rho[i] == 1) >> (gamma[i, d.vm_mypm[i]] == 0))

# Objective function
model.setObjective(quicksum(d.cluster_innerpakage_size[d.vm_mycluster[i]-1] * cal_hops(
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
    print("Result = %.4f" % model.objVal)
    print("GAP = %.4f %%" % model.MIPGap)
    print("Time = %.4f seg" % model.Runtime)

    rho_result = [1 if rho[i].X > 0.9 else 0 for i in range(d.vm_count)]
    # omiga_result = [omiga[i].X for i in range(d.vm_count)]
    gamma_result = [j for i in range(d.vm_count) for j in range(
        1, d.pm_count + 1) if gamma[i, j].X > 0.9]
    print("ρ(i): ")
    print(rho_result)
    # print("ω(i): ")
    # print(omiga_result)
    print("γ(i, j): ")
    print(gamma_result)

    solution = []
    for i in range(1, d.pm_count + 1):
        delta_b = 0
        delta_m = 0
        for j in range(d.vm_count):
            if d.vm_mypm[j] == i:
                delta_b += d.cluster_innerpakage_size[d.vm_mycluster[j]-1] * cal_hops(
                    j) * d.band_unitcost
                delta_m += d.vm_storage[j] * rho[j]

        solution.append((Uc_var[i].X, Uo[i].getValue(), Ur[i].getValue(), Ub[i].getValue(
        ), delta_b.getValue() if delta_b else 0, delta_m.getValue() if delta_m else 0))
    print("solution: ")
    for i in solution:
        print(
            f"Ucorb: ({i[0]:.4f}, {i[1]:.4f}, {i[2]:.4f}, {i[3]:.4f})\t\tdelta_bm: ({i[4]:.4f}, {i[5]:.4f})")
    print()
    for i in solution:
        print(i)

    for i in solution:
        print(calE_star(i[0:4]))
