import matplotlib.pyplot as plt
import numpy as np
from numpy import log


def euler_method(start, end, h, f, start_cond):
    n = round((end - start) / h) + 1
    x_cels = [0] * n
    y_cels = [[0] * len(f) for _ in range(n)]
    x_cels[0] = start
    y_cels[0] = start_cond
    for i in range(1, n):
        for j in range(len(f)):
            y_cels[i][j] = y_cels[i - 1][j] + h * f[j](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])
        x_cels[i] = x_cels[i - 1] + h
    return x_cels, y_cels


def improved_euler_method(start, end, h, f, start_cond):
    n = round((end - start) / h) + 1
    x_cels = [0] * n
    y_cels = [[0] * len(f) for _ in range(n)]
    x_cels[0] = start
    y_cels[0] = start_cond
    for i in range(1, n):
        for j in range(len(f)):
            y_half = [0] * len(f)
            for k in range(len(f)):
                y_half[k] = y_cels[i - 1][k] + (h / 2.0) * f[k](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])
            y_cels[i][j] = y_cels[i - 1][j] + h * f[j](x_cels[i - 1] + (h / 2.0), y_half[0], y_half[1])
        x_cels[i] = x_cels[i - 1] + h
    return x_cels, y_cels


def euler_cauchy_method(start, end, h, f, start_cond):
    n = round((end - start) / h) + 1
    x_cels = [0] * n
    y_cels = [[0] * len(f) for _ in range(n)]
    x_cels[0] = start
    y_cels[0] = start_cond
    for i in range(1, n):
        for j in range(len(f)):
            y_wave = [0] * len(f)
            for k in range(len(f)):
                y_wave[k] = y_cels[i - 1][k] + h * f[k](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])
            y_cels[i][j] = y_cels[i - 1][j] + (h / 2.0) * (
                    f[j](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1]) + f[j](x_cels[i - 1] + h, y_wave[0],
                                                                                   y_wave[1]))
        x_cels[i] = x_cels[i - 1] + h
    return x_cels, y_cels


def runge_kutta_method(start, end, h, f, start_cond):
    n = round((end - start) / h) + 1
    x_cels = [0] * n
    y_cels = [[0] * len(f) for _ in range(n)]
    x_cels[0] = start
    y_cels[0] = start_cond
    for i in range(1, n):
        k_matrix = [[0] * len(f) for _ in range(4)]
        k_matrix[0][0] = h * f[0](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])
        k_matrix[0][1] = h * f[1](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])

        k_matrix[1][0] = 2.0 * h * f[0](x_cels[i - 1] + (h / 2.0),
                                        y_cels[i - 1][0] + (k_matrix[0][0] / 2.0),
                                        y_cels[i - 1][1] + (k_matrix[0][1] / 2.0))
        k_matrix[1][1] = 2.0 * h * f[1](x_cels[i - 1] + (h / 2.0),
                                        y_cels[i - 1][0] + (k_matrix[0][0] / 2.0),
                                        y_cels[i - 1][1] + (k_matrix[0][1] / 2.0))

        k_matrix[2][0] = 2.0 * h * f[0](x_cels[i - 1] + (h / 2.0),
                                        y_cels[i - 1][0] + (k_matrix[1][0] / 4.0),
                                        y_cels[i - 1][1] + (k_matrix[1][1] / 4.0))
        k_matrix[2][1] = 2.0 * h * f[1](x_cels[i - 1] + (h / 2.0),
                                        y_cels[i - 1][0] + (k_matrix[1][0] / 4.0),
                                        y_cels[i - 1][1] + (k_matrix[1][1] / 4.0))
        k_matrix[3][0] = h * f[0](x_cels[i - 1] + h, y_cels[i - 1][0] + k_matrix[2][0] / 2.0,
                                  y_cels[i - 1][1] + k_matrix[2][1] / 2.0)
        k_matrix[3][1] = h * f[1](x_cels[i - 1] + h, y_cels[i - 1][0] + k_matrix[2][0] / 2.0,
                                  y_cels[i - 1][1] + k_matrix[2][1] / 2.0)

        y_cels[i][0] = y_cels[i - 1][0] + (sum(k_matrix[k][0] for k in range(4)) / 6.0)
        y_cels[i][1] = y_cels[i - 1][1] + (sum(k_matrix[k][1] for k in range(4)) / 6.0)
        x_cels[i] = x_cels[i - 1] + h
    return x_cels, y_cels


def adams_method(start, end, h, f, start_cond):
    n = round((end - start) / h) + 1
    x_cels = [0] * n
    y_cels = [[0] * len(f) for _ in range(n)]
    x_cels_4, y_cels_4 = runge_kutta_method(start, start + h * 3.0, h, f, start_cond)
    x_cels[0:4] = x_cels_4
    y_cels[0:4] = y_cels_4
    for i in range(4, n):
        for j in range(len(f)):
            f_4 = f[j](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])
            f_3 = f[j](x_cels[i - 2], y_cels[i - 2][0], y_cels[i - 2][1])
            f_2 = f[j](x_cels[i - 3], y_cels[i - 3][0], y_cels[i - 3][1])
            f_1 = f[j](x_cels[i - 4], y_cels[i - 4][0], y_cels[i - 4][1])
            y_cels[i][j] = y_cels[i - 1][j] + (h / 24.0) * (55.0 * f_4 - 59.0 * f_3 + 37.0 * f_2 - 9.0 * f_1)
        x_cels[i] = x_cels[i - 1] + h
    return x_cels, y_cels


def runge_romberg_method(start, end, h, function, start_cond, method, p):
    x_cels_1, y_cels_1 = method(start, end, h, function, start_cond)
    x_cels_2, y_cels_2 = method(start, end, h / 2.0, function, start_cond)
    y_cels = [0] * len(x_cels_1)
    for i in range(len(x_cels_1)):
        y_cels[i] = y_cels_2[i * 2][0] + (y_cels_2[i * 2][0] - y_cels_1[i * 1][0]) / (2.0 ** p - 1.0)
    return x_cels_1, y_cels


def main():
    start = 1.0
    end = 2.0
    y_pr = lambda x, y, z: z
    z_pr = lambda x, y, z: ((2 * x + 1) / (x * (x + 1))) * z - ((2 * x + 1) / (x ** 2 * (x + 1))) * y
    start_cond = [2.0, 4.0]
    h_po_um = 0.1
    y_resh = lambda x: (x ** 2) + x + (x * log(x))

    h = float(input() or h_po_um)
    if h <= 1e-10 or h > (end - start):
        raise ValueError(f"Step h must be between {0} and {end - start}")

    n = (end - start) / h

    if abs(n - round(n)) >= 1e-10:
        raise ValueError(f"h should completely divide the {end - start}")

    n = round(n)
    plt.figure(figsize=(10, 6))

    all_methods = {
        "Метод Эйлера": [euler_method, "red", 1.0],
        "Улучшенный метод Эйлера": [improved_euler_method, "green", 2.0],
        "Метод Коши-Эйлера": [euler_cauchy_method, "purple", 2.0],
        "Метод Рунге-Кутты четвертого порядка точности": [runge_kutta_method, "blue", 4.0],
        "Метод Адамса четвёртого порядка точности": [adams_method, "pink", 4.0]
    }

    for name, [method, color, p] in all_methods.items():
        x_res, y_res = method(start, end, h, [y_pr, z_pr], start_cond)
        x_u, y_u = runge_romberg_method(start, end, h, [y_pr, z_pr], start_cond, method, p)
        y_res_0 = []
        print(name)
        print(f"   x   |  Численное   | Аналитическое | Приближение РРР | Погрешность  |")
        print("-" * 56)
        for i in range(n + 1):
            print(
                f"{x_res[i]:.4f} | {y_res[i][0]:.10f} | {y_resh(x_res[i]):.10f}  |  {y_u[i]:.10f}   | {abs(y_resh(x_res[i]) - y_u[i]):.10f} |")
            y_res_0.append(y_res[i][0])
        print("\n\n")
        plt.plot(x_res, y_res_0, 'o-', label=name, color=color, linewidth=2)


    x = np.linspace(start, end, n)
    plt.plot(x, y_resh(x), "y--", label='Аналитическое решение $y = x^2 + x + x \ln(x)$', linewidth=2)
    plt.title("Решение дифференциального уравнения $x^2 (x + 1) y'' - x(2x+1)y' + (2x+1)y = 0$")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


if __name__ == "__main__":
    main()
