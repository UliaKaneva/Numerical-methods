from math import tan, cos, sin


def euler_method(start, end, h, f, start_cond):
    n = round((end - start) / h) + 1
    x_cels = [0] * n
    y_cels = [[0] * len(f) for _ in range(n)]
    x_cels[0] = start
    y_cels[0] = start_cond
    for i in range(1, n):
        for j in range(len(f)):
            y_cels[i][j] = y_cels[i - 1][j] + h * f[j](x_cels[i - 1], y_cels[i - 1][0], y_cels[i - 1][1])
        x_cels[i] = x_cels[i - 1] + h
    return x_cels, y_cels




def main():

    start = 0.0
    end = 1.0
    y_pr = lambda x, y, z: z
    z_pr = lambda x, y, z: -z * tan(x) - y * cos(x) ** 2
    start_cond = [0.0, 1.0]
    h_po_um = 0.1
    y_resh = lambda x: sin(sin(x))


    h = float(input() or h_po_um)
    if h < 0 or h > (end - start):
        raise ValueError(f"Step h must be between {0} and {end - start}")
    n = (end - start) / h
    if abs(n - round(n)) >= 1e-10:
        raise ValueError(f"h should completely divide the {end - start}")
    n = round(n)
    x_euler, y_euler = euler_method(start, end, h, [y_pr, z_pr], start_cond)
    print("x  y_e, y_r")
    for i in range(n + 1):
        print(f"{x_euler[i]:.3f}, {y_euler[i][-1]:.10f},")




if __name__ == "__main__":
    main()

