def read_matrix():
    n, m = map(int, input().split())
    matrix = [[float(j) for j in input().split()] for _ in range(n)]
    if any(filter(lambda x: len(x) != m, matrix)):
        raise ValueError("Количество элементов в строке матрицы не соответствует ожидаемому")
    return matrix


def read_vector():
    n = int(input())
    vector = [float(input()) for _ in range(n)]
    return vector
