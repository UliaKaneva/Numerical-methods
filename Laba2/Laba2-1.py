from math import sin, cos

Q = 0.72
A, B = 0.7, 0.80


def phi(x):
    return ((sin(x) + 0.5) / 2) ** 0.5


def f(x):
    return sin(x) - 2 * x ** 2 + 0.5


def f_fp(x):
    return cos(x) - 4 * x


def f_sp(x):
    return -sin(x) - 4


def method_simple_iteration(func, e, q, a, b):
    cof = q / (1 - q)
    x_last = func((a + b) / 2)
    count_iter = 0

    while abs(((x_current := func(x_last)) - x_last) * cof) >= e:
        count_iter += 1
        print(f"Iteration №{count_iter}: x - {x_current}; e - {x_current - x_last} ")
        x_last = x_current
    count_iter += 1
    print(f"Iteration №{count_iter}: x - {x_current}; e - {x_current - x_last} ")
    return x_current


def search_x(func, func_sp, a, b):
    x = a
    step = (b - a) / 1000
    for i in range(1000):
        if func_sp(x) * func(x) > 0:
            return x
        x += step
    raise Exception("No solution")


def method_newton(func, func_fp, func_sp, e, a, b):
    x_last = search_x(func, func_sp, a, b)
    count_iter = 0
    while abs((x_current := x_last - func(x_last) / func_fp(x_last)) - x_last) > e:
        count_iter += 1
        print(f"Iteration №{count_iter}: x - {x_current}; e - {x_current - x_last} ")
        x_last = x_current
    count_iter += 1
    print(f"Iteration №{count_iter}: x - {x_current}; e - {x_current - x_last} ")
    return x_current


def main():
    epsilon = float(input('Enter a number for epsilon: '))
    if abs(epsilon) < 1e-14 or epsilon < 0:
        raise ValueError('Epsilon cannot be <= 0')
    result = method_simple_iteration(phi, epsilon, Q, A, B)
    print(result)
    result = method_newton(f, f_fp, f_sp, epsilon, A, B)
    print(result)


if __name__ == '__main__':
    main()
