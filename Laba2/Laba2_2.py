from math import sin, cos
from par_and_func import matrix_multiplication, matrix_difference

A, B = 0.8, 0.9  # границы для x_1
C, D = 1.7, 1.8  # границы для x_2
Q = 0.975
phi = [lambda x, y: 1 + cos(y), lambda x, y: 1 + sin(x)]

resh = [0.83218792214600197821, 1.7394061793139047621]


def function(x, y):
    return [[x - cos(y) - 1],
            [y - sin(x) - 1]]


def matrix_jacoby_inv(x, y):
    det = (1 + sin(y) * cos(x))
    return [[1 / det, -sin(y) / det],
            [cos(x) / det, 1 / det]]


def error_rate(x_1, x_2):
    diff = matrix_difference(x_1, x_2)
    return max(abs(diff[0][0]), abs(diff[1][0]))





def method_simple_iteration(a, b, c, d, q, func, e):
    x_last = [(a + b) / 2, (c + d) / 2]
    x_current = [func[0](*x_last), func[1](*x_last)]
    count_iter = 0
    cof = q / (1 - q)
    pog_last = max(abs(x_current[0] - x_last[0]), abs(x_current[1] - x_last[1]))
    while pog_last * cof > e:
        count_iter += 1
        print(f"Iteration № {count_iter}: x_1 - {x_current[0]}, x_2 - {x_current[1]}; e - {pog_last}")
        x_last = x_current
        x_current = [func[0](*x_last), func[1](*x_last)]
        pog_last = max(abs(x_current[0] - x_last[0]), abs(x_current[1] - x_last[1]))
    count_iter += 1
    print(f"Iteration №{count_iter}: x_1 - {x_current[0]}, x_2 - {x_current[1]}; e - {pog_last}")
    return x_current


def method_newton(a, b, c, d, func, mt_jacoby_inv, e):
    x_last = [[(a + b) / 2],
              [(c + d) / 2]]
    x_current = matrix_difference(x_last, matrix_multiplication(mt_jacoby_inv(x_last[0][0], x_last[1][0]),
                                                                func(x_last[0][0], x_last[1][0])))
    count_iter = 0
    pog = error_rate(x_last, x_current)
    while pog > e:
        count_iter += 1
        print(f"Iteration № {count_iter}: x_1 - {x_current[0][0]}, x_2 - {x_current[1][0]}; e - {pog} ")
        x_last = x_current
        x_current = matrix_difference(x_last, matrix_multiplication(mt_jacoby_inv(x_last[0][0], x_last[1][0]),
                                                                    func(x_last[0][0], x_last[1][0])))
        pog = error_rate(x_last, x_current)
    count_iter += 1
    print(f"Iteration № {count_iter}: x_1 - {x_current[0][0]}, x_2 - {x_current[1][0]}; e - {pog} ")
    return x_current


def main():
    epsilon = float(input('Enter a number for epsilon: '))
    if abs(epsilon) < 1e-14 or epsilon < 0:
        raise ValueError('Epsilon cannot be <= 0')
    result = method_simple_iteration(A, B, C, D, Q, phi, epsilon)
    print(*result)
    print(*resh)
    result = method_newton(A, B, C, D, function, matrix_jacoby_inv, epsilon)
    print(*[i[0] for i in result])
    print(*resh)


if __name__ == "__main__":
    main()
