from math import atan as arctg
import matplotlib.pyplot as plt
import numpy as np


# Построение многочлена Лагранжа
def omega_pr(x):
    n = len(x)
    res = [1.0] * n
    for i in range(n):
        for j in range(n):
            if j == i:
                continue
            res[i] *= x[i] - x[j]
    return res


def omega(x, ind_without, x_star=None, is_string=False):
    string = ""
    count = 1.0
    if (not is_string) and (x_star is None):
        raise ValueError("Either the 'x_star' or 'is_string' parameter must be entered.")
    for i in range(len(x)):
        if i == ind_without:
            continue
        if is_string:
            if abs(x[i]) < 1e-14:
                string = string + "(x)"
            elif x[i] > 0:
                string = string + f"(x - {x[i]})"
            else:
                string = string + f"(x + {abs(x[i])})"

        if not x_star is None:
            count *= x_star - x[i]
    return count, string


def lagrange_interpolation_polynomial(x, f, is_print=False, x_star=None):
    result = 0.0
    str_res = ""
    if (not is_print) and (x_star is None):
        raise ValueError("Either the 'x_star' or 'is_print' parameter must be entered.")
    f_list = f(x)
    o_list = omega_pr(x)
    cof = [f_list[i] / o_list[i] for i in range(len(f_list))]
    for i in range(len(cof)):
        n, s = omega(x, i, x_star=x_star, is_string=is_print)
        if not x_star is None:
            result += cof[i] * n
        if is_print:
            if abs(cof[i]) < 1e-14:
                continue
            elif str_res == "":
                str_res = f"{cof[i]:.5f}{s}"
            elif cof[i] >= 0:
                str_res = f"{str_res} + {cof[i]:.5f}{s}"
            else:
                str_res = f"{str_res} - {abs(cof[i]):.5f}{s}"
    return result, str_res

# Построение многочлена Лагранжа
def count_separated_differences(x_list, f):
    sep_dif = [[0] * i for i in range(len(x_list), 0, -1)]
    for i in range(len(x_list)):
        if i == 0:
            sep_dif[0] = f(x_list)
            continue
        for j in range(len(x_list) - i):
            sep_dif[i][j] = (sep_dif[i - 1][j] - sep_dif[i - 1][j + 1]) / (x_list[j] - x_list[j + i])
    return sep_dif


def newton_interpolation_polynomial(x_list, f, x_star=None, is_print=False):
    result = 0.0
    str_res = ""
    if (not is_print) and (x_star is None):
        raise ValueError("Either the 'x_star' or 'is_print' parameter must be entered.")
    sep_dif = count_separated_differences(x_list, f)
    for i in range(len(x_list)):
        n, s = omega(x_list[:i], len(x_list), is_string=is_print, x_star=x_star)
        if not x_star is None:
            result += sep_dif[i][0] * n
        if is_print:
            if abs(sep_dif[i][0]) < 1e-14:
                continue
            elif str_res == "":
                str_res = f"{sep_dif[i][0]:.5f}{s}"
            elif sep_dif[i][0] >= 0:
                str_res = f"{str_res} + {sep_dif[i][0]:.5f}{s}"
            else:
                str_res = f"{str_res} - {abs(sep_dif[i][0]):.5f}{s}"
    return result, str_res


def main():
    x_star = -0.5
    x_a = [-3.0, -1.0, 1.0, 3.0]
    x_b = [-3.0, 0.0, 1.0, 3.0]
    f = lambda x: [arctg(x_i) for x_i in x]

    # Лагранж
    print("Интерполяционный многочлен Лагранжа")
    print("A:")
    res_num, str_res = lagrange_interpolation_polynomial(x_a, f, is_print=True, x_star=x_star)
    print(str_res)
    print(abs(res_num - arctg(x_star)))
    print("B:")
    res_num, str_res = lagrange_interpolation_polynomial(x_b, f, is_print=True, x_star=x_star)
    print(str_res)
    print(abs(res_num - arctg(x_star)))
    # Ньютон
    print("Интерполяционный многочлен Ньютон")
    print("A:")
    res_num, str_res = newton_interpolation_polynomial(x_a, f, is_print=True, x_star=x_star)
    print(str_res)
    print(abs(res_num - arctg(x_star)))
    print("B:")
    res_num, str_res = newton_interpolation_polynomial(x_b, f, is_print=True, x_star=x_star)
    print(str_res)
    print(abs(res_num - arctg(x_star)))


def fun_lag(x_a, f, x_star):
    res, r = lagrange_interpolation_polynomial(x_a, f, is_print=False, x_star=x_star)
    return res


def fun_newton(x_list, f, x_star):
    res, r = newton_interpolation_polynomial(x_list, f, is_print=False, x_star=x_star)
    return res


def graphic():
    x_star = -0.5
    x_a = [-3.0, -1.0, 1.0, 3.0]
    x_b = [-3.0, 0.0, 1.0, 3.0]
    f = lambda q: [np.arctan(x_i) for x_i in q]

    x = np.linspace(-4, 4, 400)
    y = np.arctan(x)

    # Вычисляем интерполяционные полиномы для обоих наборов точек
    p_a_lag = fun_lag(x_a, f, x)
    p_b_lag = fun_lag(x_b, f, x)
    p_a_newton = fun_newton(x_a, f, x)
    p_b_newton = fun_newton(x_b, f, x)

    # Создаем фигуру с шестью подграфиками (2 строки, 3 столбца)
    fig, axes = plt.subplots(2, 3, figsize=(12, 6))

    # Первая строка
    ax1, ax2, ax3 = axes[0]

    # Первый график для x_a (Лагранж)
    ax1.plot(x, y, label='$y = arctg(x)$', color='red', linewidth=2)
    ax1.plot(x, p_a_lag, label='$L_3(x_a)$', color='blue', linewidth=2)
    ax1.scatter(x_a, f(x_a), color='green', s=50, zorder=5, label='Узлы интерполяции')
    ax1.set_title(f'Лагранж: $x_a = {x_a}$')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linewidth=0.5)
    ax1.axvline(x=0, color='k', linewidth=0.5)

    # Второй график для x_a (Ньютон) - поменяли местами с графиком 4
    ax2.plot(x, y, label='$y = arctg(x)$', color='red', linewidth=2)
    ax2.plot(x, p_a_newton, label='$N_3(x_a)$', color='blue', linewidth=2)
    ax2.scatter(x_a, f(x_a), color='green', s=50, zorder=5, label='Узлы интерполяции')
    ax2.set_title('Ньютон: $x_a = [-3, -1, 1, 3]$')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linewidth=0.5)
    ax2.axvline(x=0, color='k', linewidth=0.5)

    # Третий график - сравнение Лагранжа и Ньютона для x_a
    ax3.plot(x, p_a_lag, label='$L_3(x_a)$', color='red', linewidth=2)
    ax3.plot(x, p_a_newton, label='$N_3(x_a)$', color='blue', linewidth=2, linestyle='--')
    ax3.scatter(x_a, f(x_a), color='green', s=50, zorder=5, label='Узлы интерполяции')
    ax3.set_title('Сравнение методов: $x_a = [-3, -1, 1, 3]$')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='k', linewidth=0.5)
    ax3.axvline(x=0, color='k', linewidth=0.5)

    # Вторая строка
    ax4, ax5, ax6 = axes[1]

    # Четвертый график для x_b (Лагранж) - поменяли местами с графиком 2
    ax4.plot(x, y, label='$y = arctg(x)$', color='red', linewidth=2)
    ax4.plot(x, p_b_lag, label='$L_3(x_b)$', color='blue', linewidth=2)
    ax4.scatter(x_b, f(x_b), color='green', s=50, zorder=5, label='Узлы интерполяции')
    ax4.set_title('Лагранж: $x_b = [-3, 0, 1, 3]$')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=0, color='k', linewidth=0.5)
    ax4.axvline(x=0, color='k', linewidth=0.5)

    # Пятый график для x_b (Ньютон)
    ax5.plot(x, y, label='$y = arctg(x)$', color='red', linewidth=2)
    ax5.plot(x, p_b_newton, label='$N_3(x_b)$', color='blue', linewidth=2)
    ax5.scatter(x_b, f(x_b), color='green', s=50, zorder=5, label='Узлы интерполяции')
    ax5.set_title('Ньютон: $x_b = [-3, 0, 1, 3]$')
    ax5.set_xlabel('x')
    ax5.set_ylabel('y')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    ax5.axhline(y=0, color='k', linewidth=0.5)
    ax5.axvline(x=0, color='k', linewidth=0.5)

    # Шестой график - сравнение Лагранжа и Ньютона для x_b
    ax6.plot(x, p_b_lag, label='$L_3(x_b)$', color='red', linewidth=2)
    ax6.plot(x, p_b_newton, label='$N_3(x_b)$', color='blue', linewidth=2, linestyle='--')
    ax6.scatter(x_b, f(x_b), color='green', s=50, zorder=5, label='Узлы интерполяции')
    ax6.set_title('Сравнение методов: $x_b = [-3, 0, 1, 3]$')
    ax6.set_xlabel('x')
    ax6.set_ylabel('y')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    ax6.axhline(y=0, color='k', linewidth=0.5)
    ax6.axvline(x=0, color='k', linewidth=0.5)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
    graphic()
