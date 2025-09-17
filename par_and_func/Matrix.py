def read_matrix():
    n, m = map(int, input("Введите 'n' и 'm': ").split())
    print("Введите матрицу")
    matrix = [[float(j) for j in input().split()] for _ in range(n)]
    if any(filter(lambda x: len(x) != m, matrix)):
        raise ValueError("Количество элементов в строке матрицы не соответствует ожидаемому")
    return matrix


def read_vector():
    n = int(input("Введите длину вектора n: "))
    print("Введите значения вектора, каждый элемент вектора на новой строке")
    vector = [float(input()) for _ in range(n)]
    return vector
