from Laba1.Laba1_2 import solving_diagonal_matrix
import matplotlib.pyplot as plt
import numpy as np


def generating_coefficient_equations(x_list, f_list):
    h_list = [x_list[i] - x_list[i - 1] for i in range(1, len(x_list))] # нахождение h_i = x_i - x_(i-1)
    n = len(x_list) - 2 # Размер системы уравнений с c_i
    matrix = [[0.0] * n for _ in range(n)]
    vector_b = [[0.0] for _ in range(n)]
    for i in range(n): # Построение трехдиагональной системы линейных уравнений, для нахождения c_i
        if i != 0: # Коэффициенты для нижнее поддиагонали
            matrix[i][i - 1] = h_list[i]
        if i != n - 1: # Коэффициенты для верхней поддиагонали
            matrix[i][i + 1] = h_list[i + 1]
        matrix[i][i] = 2 * (h_list[i] + h_list[i + 1]) # Коэффициент для главной диагонали
        # Свободный член
        vector_b[i][0] = 3 * (
                    (f_list[i + 2] - f_list[i + 1]) / h_list[i + 1] - ((f_list[i + 1] - f_list[i]) / h_list[i]))
    result = solving_diagonal_matrix(matrix, vector_b) # решение системы, метод из первой лабораторной второе задание
    # Нахождение a_i, b_i и d_i
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
    x_i = [i for i in range(len(x_list)) if x_list[i] >= x_star] # определение интервала для x_star
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

def spline_for_graphic(x_list, x_star, cof_matrix, index): # функция для построения графика сплайна
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
    f_list = [0.06334, 0.06334, 0.06334, 0.06334, 0.06334, 0.06334,
              0.06334, 0.04433, 0.03719, 0.03377,
              0.03217, 0.03172, 0.0326, 0.03501,
              0.03885, 0.04417, 0.05112, 0.05993,
              0.07087, 0.08428, 0.10056, 0.12013]

    x_list = [0.0, 0.04, 0.08, 0.16, 0.175, 0.185,
              0.19235, 0.22236, 0.23789, 0.25006,
              0.26043, 0.2697, 0.27466, 0.28836,
              0.31015, 0.34034, 0.37979, 0.42976,
              0.49183, 0.56794, 0.66027, 0.77136]
    x_star = -0.5
    # x_star = 1.5
    x_star = 0.3

    res_n, res_str = splain(x_list, f_list, x_star, is_print=True)
    print(res_n)
    print(res_str)

    # График
    x_s = []
    p_s = []
    m_p_s = []
    len_us = 0.0
    c_m = generating_coefficient_equations(x_list, f_list)
    for i in range(len(x_list) - 1):
        x_s.append(np.linspace(x_list[i] - len_us, x_list[i + 1] + len_us, 100))
        p_s.append(spline_for_graphic(x_list, x_s[i], c_m, i))


    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(x_s[0], p_s[0], label='Кубический сплайн', color='#2ca02c', linewidth=2)
    ax1.plot(x_s[0], -p_s[0], label='Кубический сплайн', color='#2ca02c', linewidth=2)
    for i in range(1, len(x_s)):
        ax1.plot(x_s[i], p_s[i], color='#2ca02c', linewidth=2)
        ax1.plot(x_s[i], -p_s[i], color='#2ca02c', linewidth=2)
    # ax1.plot(x1, p1, label='Кубический сплайн', color='#1f77b4', linewidth=2)
    # ax1.plot(x2, p2, color='#ff7f0e', linewidth=2)
    # ax1.plot(x3, p3, color='#2ca02c', linewidth=2)
    # ax1.plot(x4, p4, color='#d62728', linewidth=2)
    # ax1.scatter(x_list, f_list, color='#FF1493', s=50, zorder=5, label='Узлы интерполяции', edgecolors='black',
    #             linewidth=1)  # Ярко-розовый
    # ax1.scatter(x_star, spline_for_graphic(x_list, x_star, c_m, 7), color='#FFD700', s=80, zorder=5, label='$X^*$',
    #             edgecolors='black', linewidth=1)  # Золотой
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
