from Laba1.Laba1_1 import solve_equation
import matplotlib.pyplot as plt
import numpy as np


def approximating_polynomial(x_list, f_list, degrees=1):
    if degrees < 0:
        raise ValueError("degrees must be >= 0")
    matrix = [[0.0] * (degrees + 1) for _ in range(degrees + 1)]
    vector = [[0.0] for _ in range(degrees + 1)]
    for k in range(degrees + 1):
        for j in range(len(f_list)):
            vector[k][0] += f_list[j] * (x_list[j] ** k)
        for i in range(degrees + 1):
            for j in range(len(x_list)):
                matrix[k][i] +=  x_list[j] ** (k + i)
    a_list = [i[0] for i in solve_equation(matrix, vector)]
    return a_list

def polynomial(a_list, x):
    res = 0.0
    for i in range(len(a_list)):
        res += a_list[i] * (x ** i)
    return res

def sum_squares_errors(a_list, x_list, f_list):
    res = 0.0
    for i in range(len(f_list)):
        res += (f_list[i] - polynomial(a_list, x_list[i])) ** 2
    return res

def print_polynomial(a_list):
    res = ""
    for i in range(len(a_list)):
        if i == 0:
            s = ""
        elif i == 1:
            s = f"x"
        else:
            s = f"x^{i}"
        if res == "":
            res = f"{a_list[i]}{s}"
        elif a_list[i] > 0:
            res = f"{res} + {a_list[i]}{s}"
        else:
            res = f"{res} - {abs(a_list[i])}{s}"
    if res == 0:
        print(0)
    else:
        print(res)

def main():
    x_list = [-5.0, -3.0, -1.0, 1.0, 3.0, 5.0]
    f_list = [-1.3734, -1.249, -0.7854, 0.7854, 1.249, 1.3734]
    # x_list = [0.0, 1.7, 3.4, 5.1, 6.8, 8.5]
    # f_list = [0.0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155]
    a_list1 = approximating_polynomial(x_list, f_list)
    print_polynomial(a_list1)
    s_e = sum_squares_errors(a_list1, x_list, f_list)
    print(s_e)
    a_list2 = approximating_polynomial(x_list, f_list, degrees=2)
    print_polynomial(a_list2)
    s_e = sum_squares_errors(a_list2, x_list, f_list)
    print(s_e)
    # a_list3 = approximating_polynomial(x_list, f_list, degrees=5)
    # print_polynomial(a_list3)
    # s_e = sum_squares_errors(a_list3, x_list, f_list)
    # print(s_e)

    x = np.linspace(x_list[0] - 0.5, x_list[-1] + 0.5, 500)
    y1 = polynomial(a_list1, x)
    y2 = polynomial(a_list2, x)

    # y3 = polynomial(a_list3, x)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(x, y1, color='red', linewidth=2)
    ax1.scatter(x_list, f_list, color='green', s=40, zorder=5, label='Узлы интерполяции')
    ax1.set_title('Многочлен 1 степени')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.axvline(x=0, color='k', linewidth=0.5)

    ax2.plot(x, y2, color='red', linewidth=2)
    # ax2.plot(x, y3, color='black', linewidth=2)
    ax2.scatter(x_list, f_list, color='green', s=40, zorder=5, label='Узлы интерполяции')
    ax2.set_title('Многочлен 2 степени')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linewidth=0.5)
    ax2.axvline(x=0, color='k', linewidth=0.5)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()