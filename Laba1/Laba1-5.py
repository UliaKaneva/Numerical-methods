from par_and_func import *
import numpy as np


def norma_euclid(matrix):
    summ = sum([sum(map(lambda x: x ** 2, i)) for i in matrix])
    result = summ ** 0.5
    return result


def sign(number):
    if abs(number) < 1e-12:
        return 0
    return number / abs(number)


def qr(matrix):
    n = len(matrix)

    e_matrix = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    q_matrix = copy_matrix(e_matrix)
    a_matrix = copy_matrix(matrix)

    for i in range(n - 1):
        v_vector = [[0.0] for _ in range(n)]
        a = a_matrix[i][i]
        v_vector[i][0] = a + sign(a) * norma_euclid(take_column(a_matrix, i)[i:])
        for j in range(i + 1, n):
            v_vector[j][0] = a_matrix[j][i]
        cof = 2 / matrix_multiplication(transpose(v_vector), v_vector)[0][0]
        vv = matrix_multiplication(v_vector, transpose(v_vector))
        h_matrix = matrix_difference(e_matrix, multiplying_matrix_number(cof, vv))
        q_matrix = matrix_multiplication(q_matrix, h_matrix)
        a_matrix = matrix_multiplication(h_matrix, a_matrix)
    return q_matrix, a_matrix


def kvadr_urov(b, c):
    x1 = [0.0, 0.0]
    d = b ** 2 - 4 * c
    if d < 0:
        x1[0] = -b / 2
        x1[1] = (-d) ** 0.5 / 2
    else:
        x1[0] = (-b + d ** 0.5) / 2
        x1[1] = 0.0
    return x1


def com_check(a_past, a_current, i, epsilon):
    x_p = kvadr_urov(2 * (a_past[i][i] + a_past[i + 1][i + 1]),
                     a_past[i][i] * a_past[i + 1][i + 1] - a_past[i][i + 1] * a_past[i + 1][i])
    x_c = kvadr_urov(2 * (a_current[i][i] + a_current[i + 1][i + 1]),
                     a_current[i][i] * a_current[i + 1][i + 1] - a_current[i][i + 1] * a_current[i + 1][i])
    if x_c[1] == 0:
        return False
    diff = [[x_p[0] - x_c[0]], [x_p[1], x_p[1]]]
    if norma_euclid(diff) < epsilon:
        return True
    return False


def convergence_check(a_past, a_current, epsilon):
    n = len(a_past)
    for i in range(n - 1):
        if abs(a_current[i + 1][i]) > epsilon:
            if (i == n - 2) or (abs(a_current[i + 2][i]) < epsilon):
                if not com_check(a_past, a_current, i, epsilon):
                    return False
            return False
    return True


def qr_iteration(matrix, epsilon):
    a_past = copy_matrix(matrix)
    q_m, r_m = qr(a_past)
    a_current = matrix_multiplication(r_m, q_m)
    while not convergence_check(a_past, a_current, epsilon):
        q_m, r_m = qr(a_current)
        a_past, a_current = a_current, matrix_multiplication(r_m, q_m)
    result = []
    n = len(matrix)
    i = 0
    while i < n - 1:
        if abs(a_current[i + 1][i]) < epsilon:
            result.append([(a_current[i][i], 0.0)])
            i += 1
        else:
            x = kvadr_urov(2 * (a_current[i][i] + a_current[i + 1][i + 1]),
                           a_current[i][i] * a_current[i + 1][i + 1] - a_current[i][i + 1] * a_current[i + 1][i])
            result.append([tuple(x)])
            i += 2
    if i == n - 1:
        result.append([(a_current[i][i], 0.0)])
    return result


def main():
    matrix = read_matrix()
    n = len(matrix)
    if n != len(matrix[0]):
        raise ValueError("Матрица должна быть квадратной")
    e = float(input("Введите точность: "))

    result = qr_iteration(matrix, e)
    print("Собственные значения: ")
    for i in result:
        if abs(i[0][1]) > e:
            print(f"{i[0][0]:.6f} + i{i[0][1]:.6f}")
        else:
            print(f"{i[0][0]:.6f}")
    a = np.array(matrix)
    w, v = np.linalg.eig(a)
    print("Проверка Numpy")
    print(w)


if __name__ == "__main__":
    main()
