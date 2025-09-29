def read_matrix():
    n, m = map(int, input("Введите 'n' и 'm': ").split())
    print("Введите матрицу")
    matrix = [[float(j) for j in input().split()] for _ in range(n)]
    if any(filter(lambda x: len(x) != m, matrix)):
        raise ValueError("Количество элементов в строке матрицы не соответствует ожидаемому")
    return matrix


def print_matrix(matrix):
    for i in matrix:
        for j in i:
            print(f"{j:.4f}", end=" ")
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
