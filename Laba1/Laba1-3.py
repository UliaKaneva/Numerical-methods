from par_and_func import *


def norma_c(matrix):
    max_sum = 0
    for i in matrix:
        max_sum = max(max_sum, sum(map(abs, i)))
    return max_sum


def get_alfa_beta(matrix, vector):
    n = len(matrix)
    alfa_matrix = [i[:] for i in matrix]
    beta_vector = [i[:] for i in vector]
    for k in range(n):
        max_row = k
        max_val = abs(alfa_matrix[k][k])
        for i in range(k + 1, n):
            if abs(alfa_matrix[i][k]) > max_val:
                max_val = abs(alfa_matrix[i][k])
                max_row = i
        if max_row != k:
            alfa_matrix[k], alfa_matrix[max_row] = alfa_matrix[max_row], alfa_matrix[k]
            beta_vector[k], beta_vector[max_row] = beta_vector[max_row], beta_vector[k]
        if abs(alfa_matrix[k][k]) < 1e-12:
            raise ValueError("Матрица вырождена или близка к вырожденной")
    for i in range(n):
        a = alfa_matrix[i][i]
        beta_vector[i][0] /= a
        for j in range(n):
            if j == i:
                alfa_matrix[i][i] = 0.0
            else:
                alfa_matrix[i][j] /= -a
    return alfa_matrix, beta_vector


def simple_iteration(alfa, beta, e):
    past = [i[:] for i in beta]
    current = matrix_addition(beta, matrix_multiplication(alfa, past))
    k_iteration = 0
    while norma_c(matrix_difference(current, past)) > e:
        past, current = current, matrix_addition(beta, matrix_multiplication(alfa, current))
        k_iteration += 1
    return current, k_iteration


def seidel_method(alfa, beta, e):
    n = len(alfa)
    b_matrix = [[0] * n for _ in range(n)]
    c_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i < j:
                c_matrix[i][j] = alfa[i][j]
            else:
                b_matrix[i][j] = alfa[i][j]
    e_b_matrix = inverse_matrix(
        matrix_difference([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)], b_matrix))

    alfa_new = matrix_multiplication(e_b_matrix, c_matrix)
    beta_new = matrix_multiplication(e_b_matrix, beta)

    result, k_iteration = simple_iteration(alfa_new, beta_new, e)
    return result, k_iteration


def main():
    matrix = read_matrix()
    n = len(matrix)
    if n != len(matrix[0]):
        raise ValueError("Матрица должна быть квадратной")
    if abs(determinant(matrix)) < 1e-12:
        raise ValueError("Матрица вырожденная")
    b = read_vector()
    if len(b) != n:
        raise ValueError("Вектор должен соответствовать количеству столбцов матрицы")
    e = float(input("Введите точность:"))
    alfa, beta = get_alfa_beta(matrix, b)
    result, k_iter = simple_iteration(alfa, beta, e)
    print_matrix(result, name="Ответ с помощью метода простых итераций:")
    print(f"Количество итераций: {k_iter}")
    result_s, k_iter_s = seidel_method(alfa, beta, e)
    print_matrix(result_s, name="Ответ с помощью метода Зейделя:")
    print(f"Количество итераций: {k_iter_s}")


if __name__ == "__main__":
    main()
