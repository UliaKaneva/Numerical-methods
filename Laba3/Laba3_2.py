from Laba1.Laba1_2 import solving_diagonal_matrix
import matplotlib.pyplot as plt
import numpy as np


def generating_coefficient_equations(x_list, f_list):
    h_list = [x_list[i] - x_list[i - 1] for i in range(1, len(x_list))]
    n = len(x_list) - 2
    matrix = [[0.0] * n for _ in range(n)]
    vector_b = [[0.0] for _ in range(n)]
    for i in range(n):
        if i != 0:
            matrix[i][i - 1] = h_list[i]
        if i != n - 1:
            matrix[i][i + 1] = h_list[i + 1]
        matrix[i][i] = 2 * (h_list[i] + h_list[i + 1])
        vector_b[i][0] = 3 * (
                    (f_list[i + 2] - f_list[i + 1]) / h_list[i + 1] - ((f_list[i + 1] - f_list[i]) / h_list[i]))
    result = solving_diagonal_matrix(matrix, vector_b)
    c_list =  [0.0] + result
    a_list = [f_list[i] for i in range(len(f_list) - 1)]
    d_list = [0 for _ in range(len(f_list) - 1)]
    b_list = [0 for _ in range(len(f_list) - 1)]
    for i in range(1, len(h_list)):
        b_list[i - 1] = (f_list[i] - f_list[i - 1]) / h_list[i - 1] - h_list[i - 1] * (
                c_list[i] + 2.0 * c_list[i - 1]) / 3.0
        d_list[i - 1] = (c_list[i] - c_list[i - 1]) / (3.0 * h_list[i - 1])
    b_list[len(h_list) - 1] = (f_list[-1] - f_list[-2]) / h_list[-1] - h_list[-1] * c_list[-1] * 2.0 / 3.0
    d_list[len(h_list) - 1] = -c_list[-1] / h_list[-1] / 3.0
    cof_matrix = [a_list, b_list, c_list, d_list]
    return cof_matrix


def splain(x_list, f_list, x_star, cof_matrix=None, is_print=False):
    result = 0.0
    res_str = ""
    if cof_matrix is None:
        cof_matrix = generating_coefficient_equations(x_list, f_list)
    x_i = [i for i in range(len(x_list)) if x_list[i] >= x_star]
    if len(x_i) == 0 or len(x_i) == len(x_list):
        return result, res_str
    index = x_i[0] - 1
    if is_print:
        s = ""
        if abs(x_list[index]) < 1e-14:
            s = "x"
        elif x_list[index] > 0:
            s = f"(x - {x_list[index]})"
        else:
            s = f"(x + {abs(x_list[index])})"
        for i in range(4):
            if abs(cof_matrix[i][index]) < 1e-14:
                continue
            if i == 0:
                res_str = f"{cof_matrix[i][index]:.5f}"
            elif cof_matrix[i][index] > 0:
                if i == 1:
                    res_str = f"{res_str} + {cof_matrix[i][index]:.5f}{s}"
                else:
                    res_str = f"{res_str} + {cof_matrix[i][index]:.5f}{s}^{i}"
            else:
                if i == 1:
                    res_str = f"{res_str} - {abs(cof_matrix[i][index]):.5f}{s}"
                else:
                    res_str = f"{res_str} - {abs(cof_matrix[i][index]):.5f}{s}^{i}"
    n = x_star - x_list[index]
    for i in range(4):
        result += cof_matrix[i][index] * (n ** i)

    return result, res_str

def spline_for_graphic(x_list, x_star, cof_matrix, index):
    result = 0.0
    n = x_star - x_list[index]
    for i in range(4):
        result += cof_matrix[i][index] * (n ** i)
    return result

def main():
    x_list = [-3.0, -1.0, 1.0, 3.0, 5.0]
    f_list = [-1.2490, -0.78540, 0.78540, 1.2490, 1.3734]
    x_star = -0.5
    res_n, res_str = splain(x_list, f_list, x_star, is_print=True)
    print(res_n)
    print(res_str)

    x_s = []
    p_s = []
    len_us = 0.5
    c_m = generating_coefficient_equations(x_list, f_list)
    for i in range(len(x_list) - 1):
        x_s.append(np.linspace(x_list[i] - len_us, x_list[i + 1] + len_us, 100))
        p_s.append(spline_for_graphic(x_list, x_s[i], c_m, i))


    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(x_s[0], p_s[0], label='Кубический сплайн', color='#1f77b4', linewidth=2)
    ax1.plot(x_s[1], p_s[1], color='#ff7f0e', linewidth=2)
    ax1.plot(x_s[2], p_s[2], color='#2ca02c', linewidth=2)
    ax1.plot(x_s[3], p_s[3], color='#d62728', linewidth=2)
    ax1.scatter(x_list, f_list, color='#FF1493', s=50, zorder=5, label='Узлы интерполяции', edgecolors='black',
                linewidth=1)
    ax1.set_title(f'Сплайн')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.axvline(x=0, color='k', linewidth=0.5)
    plt.axis('equal')

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
