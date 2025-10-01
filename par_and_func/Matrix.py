def read_matrix():
    n, m = map(int, input("Введите 'n' и 'm': ").split())
    print("Введите матрицу")
    matrix = [[float(j) for j in input().split()] for _ in range(n)]
    if any(filter(lambda x: len(x) != m, matrix)):
        raise ValueError("Количество элементов в строке матрицы не соответствует ожидаемому")
    return matrix


def print_matrix(matrix, name=""):
    if name != "":
        print(name)
    for i in matrix:
        for j in i:
            print(f"{j:.6f}", end=" ")
        print()


def read_vector():
    n = int(input("Введите длину вектора n: "))
    print("Введите значения вектора, каждый элемент вектора на новой строке")
    vector = [[float(input())] for _ in range(n)]
    return vector


def transpose(matrix):
    n, m = len(matrix), len(matrix[0])
    result = [[matrix[j][i] for j in range(n)] for i in range(m)]
    return result


def take_column(matrix, index):
    result = [[i[index]] for i in matrix]
    return result


def copy_matrix(matrix):
    return [i[:] for i in matrix]


def matrix_multiplication(matrix_1, matrix_2):
    if len(matrix_1[0]) != len(matrix_2):
        raise ValueError("Количество столбцов в первой матрице должно быть равно количеству строк второй матрице!")
    result = [[0 for _ in matrix_2[0]] for _ in matrix_1]
    for i in range(len(matrix_1)):
        for j in range(len(matrix_2[0])):
            for k in range(len(matrix_2)):
                result[i][j] += matrix_1[i][k] * matrix_2[k][j]
            if abs(result[i][j]) < 1e-12:
                result[i][j] = 0.0
    return result


def multiplying_matrix_number(number, matrix):
    result = [[matrix[i][j] * number for j in range(len(matrix[0]))] for i in range(len(matrix))]
    return result


def matrix_addition(matrix_1, matrix_2):
    if len(matrix_1) != len(matrix_2) or len(matrix_1[0]) != len(matrix_2[0]):
        raise ValueError("Матрицы должны быть одинакового размера!")
    result = [[matrix_1[i][j] + matrix_2[i][j] for j in range(len(matrix_1[0]))] for i in range(len(matrix_1))]
    return result


def matrix_difference(matrix_1, matrix_2):
    result = matrix_addition(matrix_1, multiplying_matrix_number(-1, matrix_2))
    return result


def lup(matrix):
    n = len(matrix)
    if n <= 0:
        raise ValueError("Матрица вырождена")

    count_swaps = 0
    a_matrix = [row[:] for row in matrix]
    l_matrix = [[0.0] * n for _ in range(n)]
    u_matrix = [[0.0] * n for _ in range(n)]
    p_matrix = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for k in range(n):
        max_row = k
        max_val = abs(a_matrix[k][k])
        for i in range(k + 1, n):
            if abs(a_matrix[i][k]) > max_val:
                max_val = abs(a_matrix[i][k])
                max_row = i

        if max_row != k:
            count_swaps += 1
            a_matrix[k], a_matrix[max_row] = a_matrix[max_row], a_matrix[k]
            p_matrix[k], p_matrix[max_row] = p_matrix[max_row], p_matrix[k]
            for j in range(k):
                l_matrix[k][j], l_matrix[max_row][j] = l_matrix[max_row][j], l_matrix[k][j]

        if abs(a_matrix[k][k]) < 1e-12:
            raise ValueError("Матрица вырождена или близка к вырожденной")

        l_matrix[k][k] = 1.0

        for j in range(k, n):
            u_matrix[k][j] = a_matrix[k][j]
            for m in range(k):
                u_matrix[k][j] -= l_matrix[k][m] * u_matrix[m][j]

        for i in range(k + 1, n):
            l_matrix[i][k] = a_matrix[i][k]
            for m in range(k):
                l_matrix[i][k] -= l_matrix[i][m] * u_matrix[m][k]
            l_matrix[i][k] /= u_matrix[k][k]

    return l_matrix, u_matrix, p_matrix, count_swaps


def determinant(matrix):
    l_matrix, u_matrix, p_matrix, count_swaps = lup(matrix)
    result = -1 ** (count_swaps % 2)
    for i in range(len(matrix)):
        result *= u_matrix[i][i]
    return result


def inverse_matrix(matrix):
    det = determinant(matrix)
    if abs(det) < 1e-12:
        raise ValueError("Матрица вырожденная")
    n = len(matrix)
    l_matrix, u_matrix, p_matrix, count_swap = lup(matrix)
    x = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(i):
                x[i][k] -= x[j][k] * l_matrix[i][j]
            if abs(x[i][k]) < 1e-10:
                x[i][k] = 0.0
    for k in range(n):
        for i in range(n - 1, -1, -1):
            for j in range(n - 1, i, -1):
                x[i][k] -= x[j][k] * u_matrix[i][j]
            x[i][k] /= u_matrix[i][i]
            if abs(x[i][k]) < 1e-10:
                x[i][k] = 0.0
    return matrix_multiplication(x, p_matrix)
