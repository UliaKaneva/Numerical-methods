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
    # x_list = [0.0, 1.0, 2.0, 3.0, 4.0]
    # f_list = [0.0, 1.8415, 2.9093, 3.1411, 3.2432]
    x_star = -0.5
    # x_star = 1.5

    res_n, res_str = splain(x_list, f_list, x_star, is_print=True)
    print(res_n)
    print(res_str)

    # График
    x1 = np.linspace(x_list[0], x_list[1], 100)
    x2 = np.linspace(x_list[1], x_list[2], 100)
    x3 = np.linspace(x_list[2], x_list[3], 100)
    x4 = np.linspace(x_list[3], x_list[4], 100)
    c_m = generating_coefficient_equations(x_list, f_list)
    p1 = spline_for_graphic(x_list, x1, c_m, 0)
    p2 = spline_for_graphic(x_list, x2, c_m, 1)
    p3 = spline_for_graphic(x_list, x3, c_m, 2)
    p4 = spline_for_graphic(x_list, x4, c_m, 3)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(x1, p1, label='Кубический сплайн', color='red', linewidth=2)
    ax1.plot(x2, p2, color='red', linewidth=2)
    ax1.plot(x3, p3, color='red', linewidth=2)
    ax1.plot(x4, p4, color='red', linewidth=2)
    ax1.scatter(x_list, f_list, color='green', s=40, zorder=5, label='Узлы интерполяции')
    ax1.scatter(x_star, spline_for_graphic(x_list, x_star, c_m, 1), color='blue', s=50, zorder=5, label='$X^*$')
    ax1.set_title(f'Сплайн')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.axvline(x=0, color='k', linewidth=0.5)

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
