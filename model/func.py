from sympy import diff, symbols, log
import numpy as np
import time
import random
import math
random.seed(42)


# P是0或者1,布尔值


def run(p_f, p_l, m_b, m_c, mi_lk, P, m, tau, x, cx3, cx4):
    # 定义开始时间
    time_start = time.time()

    # 定义变量
    x1, x2, x3, x4, t = symbols('x1 x2 x3 x4 t')
    # ((p_f - p_l) * 1000 * x2 * x1 + p_l + x4 * m_b + x3 * P)
    # + x2 * m_c * (x3 / mi_lk) * P
    # 定义目标函数
    func = t * ((p_f - p_l) * 1000 * x2 * x1 + p_l + x3 * cx3 + x4 * cx4) - log(
        -(tau - (1 / (6 * m) * (3 + log(1 - 1 / 2 * x1) +
                                log(1 - 1 / 2 * x2) + log(1 - 1 / 2 * x3)) + 1 / (2 * m) * (
            1 + log(1 - 1 / 2 * x4))))) - log(x1) - log(1 - x1) - log(x2) - log(1 - x2) - log(x3) - \
        log(1 - x3) - log(x4) - log(1 - x4)

    # 求导
    de_x1 = diff(func, x1, 1)  # 对x1求一阶导
    de_x2 = diff(func, x2, 1)  # 对x2求一阶导
    de_x3 = diff(func, x3, 1)  # 对x3求一阶导
    de_x4 = diff(func, x4, 1)  # 对x4求一阶导

    de_x1_x1 = diff(func, x1, x1)
    de_x1_x2 = diff(func, x1, x2)
    de_x1_x3 = diff(func, x1, x3)
    de_x1_x4 = diff(func, x1, x4)

    de_x2_x1 = diff(func, x2, x1)
    de_x2_x2 = diff(func, x2, x2)
    de_x2_x3 = diff(func, x2, x3)
    de_x2_x4 = diff(func, x2, x4)

    de_x3_x1 = diff(func, x3, x1)
    de_x3_x2 = diff(func, x3, x2)
    de_x3_x3 = diff(func, x3, x3)
    de_x3_x4 = diff(func, x3, x4)

    de_x4_x1 = diff(func, x4, x1)
    de_x4_x2 = diff(func, x4, x2)
    de_x4_x3 = diff(func, x4, x3)
    de_x4_x4 = diff(func, x4, x4)

    # 计算目标函数在x处的一阶导数
    def gradient(x_1, x_2, x_3, x_4, _t):
        # 可以选择使用实部或虚部
        j1 = float(de_x1.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        j2 = float(de_x2.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        j3 = float(de_x3.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        j4 = float(de_x4.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        return np.asmatrix([j1, j2, j3, j4]).T

    def generate_vector(a):
        for i in range(len(a)):
            random_factor = random.uniform(0.5, 1.8)
            a[i, 0] *= random_factor
        while any(element[0] > 1.0 for element in a):
            a = generate_vector(a)
        return a

    # 求二阶导数矩阵(海塞矩阵)的逆矩阵
    def inv_hessian(x_1, x_2, x_3, x_4, _t):
        x_1, x_2, x_3, x_4 = float(x_1), float(x_2), float(x_3), float(x_4)
        h11 = float(de_x1_x1.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h12 = float(de_x1_x2.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h13 = float(de_x1_x3.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h14 = float(de_x1_x4.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])

        h21 = float(de_x2_x1.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h22 = float(de_x2_x2.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h23 = float(de_x2_x3.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h24 = float(de_x2_x4.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])

        h31 = float(de_x3_x1.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h32 = float(de_x3_x2.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h33 = float(de_x3_x3.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h34 = float(de_x3_x4.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])

        h41 = float(de_x4_x1.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h42 = float(de_x4_x2.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h43 = float(de_x4_x3.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])
        h44 = float(de_x4_x4.subs(
            {x1: x_1, x2: x_2, x3: x_3, x4: x_4, t: _t}).as_real_imag()[0])

        H = np.asmatrix([[h11, h12, h13, h14], [h21, h22, h23, h24], [h31, h32, h33, h34],
                         [h41, h42, h43, h44]])
        epsilon = 1e-6  # 可以根据实际情况调整这个值
        # 如果矩阵不可逆，添加正则化
        while np.linalg.cond(H) > 1 / epsilon:
            H = H + epsilon * np.eye(H.shape[0])
        H_1 = np.linalg.inv(H)
        return H_1

    t_value = 0.01
    # x是牛顿法的初始迭代值,这样表示的是列向量。[10,10]是行向量
    # x = np.array([[0.7], [0.1], [0.1], [0.1]])
    x = np.array([[i] for i in x])
    t_iter_cnt = 0  # 记录迭代的次数
    eps = 0.00001  # 迭代停止的误差
    consume_time = 0
    while t_iter_cnt < 1:
        t_iter_cnt += 1
        # 明确目标函数
        # 迭代停止的控制条件
        iter_cnt = 0
        while iter_cnt < 7:
            iter_cnt += 1
            # print(iter_cnt)
            J = gradient(x[0, 0], x[1, 0], x[2, 0], x[3, 0], t_value)
            H_1 = inv_hessian(x[0, 0], x[1, 0], x[2, 0], x[3, 0], t_value)
            direction = -np.matmul(H_1, J)
            step_size = 0.2
            # 更新迭代点
            x_new = x + step_size * direction

            error = np.linalg.norm(J)  # 求二范数,判断迭代效果
            x = x_new
            if error < eps:
                break
        if t_value <= eps:
            break
        t_value = t_value / 10

    time_end = time.time()
    consume_time = time_end - time_start

    # 生成向量
    result_vector = generate_vector(x)
    # 保留每个元素小数点后两位
    rounded_result = np.round(result_vector, decimals=2)
    # * object_function_value = (p_f - p_l) * x[1, 0] * x[0, 0] + p_l + x[3, 0] * m_b + x[1, 0] * m_c * (
    # * x[2, 0] / mi_lk) * P
    object_function_value = None

    # print(rounded_result)
    return rounded_result, object_function_value, consume_time


def calRmj(u_m_j):
    """
    期待为0-1之间的数值
    :param u_m_j:
    :return:
    """
    return 1+math.log(1-(1-1/math.e)*u_m_j)


def getE(uc, uo, ur, ub):
    tmplist = [uc, uo, ur, ub]
    return min([calRmj(i) for i in tmplist])


def _calE_star_Ucorb(u_corb, cx3, cx4):
    # p_f, p_l, m_b, m_c, mi_lk, P, m, tau
    # p_f = 135  #满负载能耗，random(113，135)
    # p_l = 58.4 #空负载能耗，random(42.3,93.7)
    # m_b = 4    #当前pm到其他15个pm跳数的平均值。
    # m_c = 20000 #cpu
    # mi_lk = 1  #链路带宽，到其他pm的最小带宽值
    # rou = 0/1  #0，1算取val的最小的值
    # m = 16     #16个虚拟机
    # tau = 0.2  #
    # u_corb          #当时的pm的负载率
    vec = []
    vec.append(run(125, 60, None, None, None, None, 7, 0.2, u_corb, cx3, cx4))
    # * vec.append(run(random.randint(113, 135), random.randint(
    # * 423, 937)/10,  4, 20000, 1, 1, 16, 0.3, u_corb))  # 返回U和目标函数值和运行时间
    # * vec.append(run(random.randint(113, 135), random.randint(
    # * 423, 937)/10,  4, 20000, 1, 0, 16, 0.3, u_corb))  # 返回U和目标函数值和运行时间
    # * vec = sorted(vec, key=lambda x: x[1])
    retval = None
    for ans in vec:
        res = [i[0] for i in ans[0]]
        flag = True
        for i in res:
            if i < 1e-5 or i > 1:
                flag = False
        if flag == True:
            retval = res
            break

    return getE(retval[0], retval[1], retval[2], retval[3])


def calE_star(u_corb, cx3, cx4):
    # cx3: x3的系数
    if min(u_corb) < 0.00001:
        return 1
    if max(u_corb) > 0.99999:
        return 0
    # print("Now is cal e star")
    retval = None
    while retval == None:
        try:
            retval = _calE_star_Ucorb(u_corb, cx3, cx4)
        except TypeError as e:
            pass
    return retval


if __name__ == '__main__':
    print(calE_star([0.1, 0.9, 0.9, 0.8]))
