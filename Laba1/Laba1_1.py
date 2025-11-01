from par_and_func import matrix_multiplication, read_matrix, read_vector, transpose, print_matrix


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


def solve_equation(matrix, vector):
    n = len(matrix)
    if (n != len(matrix[0])) or (n != len(vector)):
        raise ValueError("Матрица должна быть квадратной и соответствовать количеству не известных")
    l_matrix, u_matrix, p_matrix, count_swaps = lup(matrix)
    b = matrix_multiplication(p_matrix, vector)

    # Ly=b
    y = [i[:] for i in b]
    for i in range(n):
        for j in range(i):
            y[i][0] -= y[j][0] * l_matrix[i][j]

    # Ux = y
    x = [i[:] for i in y]
    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            x[i][0] -= x[j][0] * u_matrix[i][j]
        x[i][0] /= u_matrix[i][i]
    return x


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


def main():
    matrix = read_matrix()
    n = len(matrix)
    if len(matrix) != len(matrix[0]):
        raise ValueError("Матрица должна быть квадратной")
    b = read_vector()
    if len(b) != n:
        raise ValueError("Вектор должен соответствовать количеству столбцов матрицы")

    l_matrix, u_matrix, p_matrix, count_swap = lup(matrix)
    print("L")
    print_matrix(l_matrix)
    print("U")
    print_matrix(u_matrix)
    print("p")
    print_matrix(p_matrix)

    result = transpose(solve_equation(matrix, b))
    print("Ответ: ", end="")
    for i in range(len(result[0])):
        print(f"x{i + 1} {result[0][i]:.4f}", end=" ")
    print()
    print(f"Определитель: {determinant(matrix):.4f}")
    inv_matrix = inverse_matrix(matrix)
    print_matrix(inv_matrix, name="Обратная матрица")
    test = matrix_multiplication(matrix, inv_matrix)
    print_matrix(test, name="Проверка:")


if __name__ == "__main__":
    main()
