def first_derivative(x_list, f_list, x_star):
    index = -1
    for i in range(len(x_list) - 1):
        if x_list[i] <= x_star <= x_list[i + 1]:
            index = i
            break
    if index == -1:
        raise ValueError(f"{x_star} is not in [{x_list[0]}, {x_list[-1]}]")
    if index == len(x_list) - 2 or index == 0:
        return (f_list[index + 1] - f_list[index]) / (x_list[index + 1] - x_list[index])
    b = (f_list[index + 1] - f_list[index]) / (x_list[index + 1] - x_list[index])
    k1 = (f_list[index + 2] - f_list[index + 1]) / (x_list[index + 2] - x_list[index + 1])
    k2 = (f_list[index + 1] - f_list[index]) / (x_list[index + 1] - x_list[index])
    k = (k1 - k2) / (x_list[index + 2] - x_list[index])
    return b + k * (2 * x_star - x_list[index] - x_list[index + 1])


def second_derivative(x_list, f_list, x_star):
    index = -1
    for i in range(len(x_list) - 1):
        if x_list[i] < x_star <= x_list[i + 1]:
            index = i
            break

    if index == -1:
        raise ValueError(f"{x_star} is not in [{x_list[0]}, {x_list[-1]}]")
    ans1 = (f_list[index + 2] - f_list[index + 1]) / (x_list[index + 2] - x_list[index + 1])
    ans2 = (f_list[index + 1] - f_list[index]) / (x_list[index + 1] - x_list[index])
    return 2 * (ans1 - ans2) / (x_list[index + 2] - x_list[index])


def main():
    x_list = [0.0, 0.5, 1.0, 1.5, 2.0]
    f_list = [0.0, 0.97943, 1.8415, 2.4975, 2.9093]
    x_star = 1.0
    # x_list = [0.0, 0.1, 0.2, 0.3, 0.4]
    # f_list = [1.0, 1.1052, 1.2214, 1.3499, 1.4918]
    # x_star = 0.2
    # x_list = [-1.0, 0.0, 1.0, 2.0, 3.0]
    # f_list = [-1.7854, 0.0, 1.7854, 3.1071, 4.249]
    # x_star = 1.0

    f_d = first_derivative(x_list, f_list, x_star)
    print(f"Первая производная в точке {x_star}: {f_d:.5f}")
    s_d = second_derivative(x_list, f_list, x_star)
    print(f"Вторая производная в точке {x_star}: {s_d:.5f}")


if __name__ == "__main__":
    main()
