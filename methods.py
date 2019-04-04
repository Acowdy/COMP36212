import math
import numpy

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
    y = numpy.zeros(n+1)
    x = numpy.zeros(n+1)
    x[0] = x0
    x[1] = x[0] + h
    y[0] = y0
    y[1] = forward_euler(f, x[0], x[1], y[0], 1)
    for n in range(1, n):
        f0 = f(x[n], y[n])
        f1 = f(x[n+1], y[n+1])
        y[n+1] = y[n] + (h / 2) * (3 * f1 - f0)
        x[n+1] = x[n]
    return y


if __name__ == "__main__":
    print(forward_euler(lambda x, y: y, 0, 4, 1, 1))
    print(ab2(lambda x, y: y, 0, 2, 1, 0.5)[-1])
