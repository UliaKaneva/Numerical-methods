from Laba4_1 import runge_kutta_method, runge_romberg_method
from numpy import exp
import matplotlib.pyplot as plt
import numpy as np
from Laba1.Laba1_1 import solve_equation


def phi(start, end, h, f, start_cond, n):
    x_res, y_res = runge_kutta_method(start, end, h, f, [n, start_cond])
    return y_res[-1][0] - y_res[-1][1] / 2.0


def shooting_method(start, end, h, f, params):
    start_cond = params[0]
    list_parm = params[1:3]
    error = params[3]
    phi_list = [phi(start, end, h, f, start_cond, i) for i in list_parm]
    while abs(phi_list[-1]) > error:
        n = list_parm[-1] - ((list_parm[-1] - list_parm[-2]) / (phi_list[-1] - phi_list[-2])) * phi_list[-1]
        list_parm.append(n)
        phi_list.append(phi(start, end, h, f, start_cond, n))
    x_res, y_res = runge_kutta_method(start, end, h, f, [list_parm[-1], start_cond])
    return x_res, y_res


def finite_difference_method(start, end, h, f, params):
    n = round((end - start) / h) + 1
    p, q = params[0], params[1]
    x_cels = [start + h * i for i in range(n)]
    vector = [[0.0] for _ in range(n)]
    matrix = [[0.0] * n for _ in range(n)]

    matrix[0][0] = -3.0
    matrix[0][1] = 4.0
    matrix[0][2] = -1.0
    vector[0][0] = 6.0 * exp(1.0) * h

    matrix[n - 1][n - 3] = 1.0
    matrix[n - 1][n - 2] = -4.0
    matrix[n - 1][n - 1] = 3.0 - 4.0 * h

    for i in range(1, n - 1):
        matrix[i][i - 1] = 1.0 - ((p(x_cels[i]) * h) / 2.0)
        matrix[i][i + 1] = 1.0 + ((p(x_cels[i]) * h) / 2.0)
        matrix[i][i] = -2.0 + (h ** 2) * q(x_cels[i])

    y_res = solve_equation(matrix, vector)

    return x_cels, y_res


def main():
    start = 1.0
    end = 2.0
    y_pr = lambda x, y, z: z
    z_pr = lambda x, y, z: ((2 * x + 1) / x) * z - ((x + 1) / x) * y
    h_po_um = 0.1
    y_resh = lambda x: exp(x) * x ** 2

    h = float(input() or h_po_um)
    if h <= 1e-10 or h > (end - start):
        raise ValueError(f"Step h must be between {0} and {end - start}")

    n = (end - start) / h

    if abs(n - round(n)) >= 1e-10:
        raise ValueError(f"h should completely divide the {end - start}")

    n = round(n)

    plt.figure(figsize=(10, 6))

    all_methods = {
        "Метод стрельбы": [shooting_method, "red", 4.0, [3.0 * exp(1.0), 1.0, 2.0, 0.001]],
        "Конечно-разностный метод": [finite_difference_method, "green", 2.0,
                                     [lambda x: -(2 * x + 1) / x, lambda x: (x + 1) / x]],
    }

    for name, [method, color, p, start_cond] in all_methods.items():
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
    plt.plot(x, y_resh(x), "y--", label='Аналитическое решение $y = x^2e^x$', linewidth=2)
    plt.title("Решение дифференциального уравнения $xy'' - x(2x+1)y' + (x+1)y = 0$")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.show()


if __name__ == "__main__":
    main()
