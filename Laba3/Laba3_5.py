def rectangle_method(start, end, step, function):  # Порядок точности p = 2
    res = 0.0
    last = start
    current = start + step
    while current <= end:
        res += function((current + last) / 2) * step
        last = current
        current = current + step
    return res


def trapezoid_method(start, end, step, function):  # Порядок точности p = 2
    res = 0.0
    last = start
    current = start + step
    while current <= end:
        res += (function(current) + function(last)) * step
        last = current
        current = current + step
    return res / 2


def simpson_method(start, end, step, function):  # Порядок точности p = 4
    res = 0.0
    last = start
    current = start + 2 * step
    while current <= end:
        res += (function(last) + 4 * function(last + step) + function(current)) * step
        last = current
        current = current + 2 * step
    return res / 3


def runge_romberg_method(start, end, step1, step2, function, method, p):
    h1, h2 = step1, step2
    if step1 > step2:
        h1, h2 = step2, step1
    k = h2 / h1
    y1 = method(start, end, h1, function)
    y2 = method(start, end, h2, function)
    delta = (y2 - y1) / (1 - k ** p)
    return abs(delta)


def main():
    a, b = 0, 2
    h1 = 0.5
    h2 = 0.25
    f = lambda x: (x ** 2) / (x ** 2 + 16)
    # a, b = -1, 1
    # h1 = 0.5
    # h2 = 0.25
    # f = lambda x: x / ((3 * x + 4) ** 2)
    resh = 0.1454095639967755
    print("Метод прямоугольников")
    print(f"С шагом h = {h1}: {rectangle_method(a, b, h1, f):.8f}")
    print(f"С шагом h = {h2}: {rectangle_method(a, b, h2, f):.8f}")
    print(f"Погрешность для h = {h2}: {runge_romberg_method(a, b, h1, h2, f, rectangle_method, 2):.8f}")
    print(f"Разница с решением: решение - {resh:.8f}, h = {h2}: {abs(rectangle_method(a, b, h2, f) - resh):.8f}",
          end="\n\n")

    print("Метод трапеций")
    print(f"С шагом h = {h1}: {trapezoid_method(a, b, h1, f):.8f}")
    print(f"С шагом h = {h2}: {trapezoid_method(a, b, h2, f):.8f}")
    print(f"Погрешность для h = {h2}: {runge_romberg_method(a, b, h1, h2, f, trapezoid_method, 2):.8f}")
    print(f"Разница с решением: решение - {resh:.8f}, h = {h2}: {abs(trapezoid_method(a, b, h2, f) - resh):.8f}",
          end="\n\n")

    print("Метод Симпсона")
    print(f"С шагом h = {h1}: {simpson_method(a, b, h1, f):.8f}")
    print(f"С шагом h = {h2}: {simpson_method(a, b, h2, f):.8f}")
    print(f"Погрешность для h = {h2}: {runge_romberg_method(a, b, h1, h2, f, simpson_method, 4):.8f}")
    print(f"Разница с решением: решение - {resh:.8f}, h = {h2}: {abs(simpson_method(a, b, h2, f) - resh):.8f}",
          end="\n\n")


if __name__ == '__main__':
    main()
