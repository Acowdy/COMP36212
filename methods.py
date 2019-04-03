import math

def forward_euler(f, x0, xn, y0, h):
    """
    Solve an IVP y' = f(x, y) using the forward Euler method where
    x0 is the initial value of x, xn is the final value of x, y0 is the
    initial value of y and h is the step size.
    """
    n = math.floor((xn - x0) / h)
    x = x0
    y = y0
    for _ in range(n):
        y += h * f(x, y)
        x += h
    return y


def ab2(f, x0, xn, y0, h):
    """
    Solve an IVP y' = f(x, y) using the AB2 method where
    x0 is the initial value of x, xn is the final value of x, y0 is the
    initial value of y and h is the step size.
    """
    n = math.floor((xn - x0) / h)
    x1 = x0 + h
    y1 = forward_euler(f, x0, x1, y0, 1)
    for _ in range(n-1):
        f0 = f(x0, y0)
        f1 = f(x1, y1)
        y0 = y1
        y1 += (h / 2) * (3 * f1 - f0)
        x0 = x1
        x1 += h
    return y1


if __name__ == "__main__":
    print(forward_euler(lambda x, y: y, 0, 4, 1, 1))
    print(ab2(lambda x, y: y, 0, 2, 1, 0.5))
