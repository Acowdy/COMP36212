from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np

Result = namedtuple("Result", ["S", "E", "ES", "P", "kf", "kr", "kc", "h", "T"])

def forward_euler(S0, E0, ES0, P0, kf, kr, kc, h, T):
    total_steps = int(T/h + 1)
    S = np.zeros(total_steps)
    E = np.zeros(total_steps)
    ES= np.zeros(total_steps)
    P = np.zeros(total_steps)
    S[0] = S0
    E[0] = E0
    ES[0]= ES0
    P[0] = P0

    t = 0
    i = 0

    dS = lambda i : -kf*E[i]*S[i] + kr*ES[i]
    dE = lambda i : -kf*E[i]*S[i] + kr*ES[i] + kc*ES[i]
    dES= lambda i :  kf*E[i]*S[i] - kr*ES[i] - kc*ES[i]
    dP = lambda i :                            kc*ES[i]

    while i+1 < total_steps:
        S[i+1] = S[i] + h*dS(i)
        E[i+1] = E[i] + h*dE(i)
        ES[i+1]= ES[i]+ h*dES(i)
        P[i+1] = P[i] + h*dP(i)
        t += h
        i += 1

    return Result(S, E, ES, P, kf, kr, kc, h, T)

def ab2_method(S0, E0, ES0, P0, kf, kr, kc, h, T):
    total_steps = int(T/h + 1)
    S = np.zeros(total_steps)
    E = np.zeros(total_steps)
    ES= np.zeros(total_steps)
    P = np.zeros(total_steps)
    S[0] = S0
    E[0] = E0
    ES[0]= ES0
    P[0] = P0

    dS = lambda i : -kf*E[i]*S[i] + kr*ES[i]
    dE = lambda i : -kf*E[i]*S[i] + kr*ES[i] + kc*ES[i]
    dES= lambda i :  kf*E[i]*S[i] - kr*ES[i] - kc*ES[i]
    dP = lambda i :                            kc*ES[i]


    # First we need to find the first value for the concentrations
    # by doing forward euler at a step size of <h^2 to the first step
    fe_result = forward_euler(S0, E0, ES0, P0, kf, kr, kc, h**2, h)

    S[1] = fe_result.S[-1]
    E[1] = fe_result.E[-1]
    ES[1]= fe_result.ES[-1]
    P[1] = fe_result.P[-1]

    t = 0
    i = 0

    while i+2 < total_steps:
        S [i+2] = (h/2)*(3* dS(i+1)-dS (i)) +  S[i+1]
        E [i+2] = (h/2)*(3* dE(i+1)-dE (i)) +  E[i+1]
        ES[i+2] = (h/2)*(3*dES(i+1)-dES(i)) + ES[i+1]
        P [i+2] = (h/2)*(3* dP(i+1)-dP (i)) +  P[i+1]
        t += h
        i += 1

    return Result(S, E, ES, P, kf, kr, kc, h, T)

def plot_result(res):
    x_labels = [x*res.h for x in range(int(res.T/res.h + 1))]
    plt.plot(x_labels, res.S)
    plt.plot(x_labels, res.E)
    plt.plot(x_labels, res.ES)
    plt.plot(x_labels, res.P)
    plt.show()


result = ab2_method(500, 200, 0, 0, 0.1, 0.05, 0.01, 0.01, 500)

plot_result(result)
