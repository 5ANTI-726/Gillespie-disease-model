from random import uniform
from operator import add
from math import log
import matplotlib.pyplot as plt


def uni():
    #Input: None
    #Output:
    return uniform(0,1)

def V_population(S, V, t):
    #Input:
    #Output:
    V_start = V[-1]
    case = 0

    #Percentage of final vaccinated people
    vax = 0.3
    #Start date of vaccination campaing
    start = 360
    #Duration of vaccination campaign
    duration = 365

    if t[-1] < start:
        V.append(0)
        case = 1
    elif t[-1] >= start and t[-1] < (start + duration):
        t0 = t[-1] - start
        expected_V = N*vax*t0/duration
        new_V = int(expected_V-V[-1])

        if S[-1] > new_V:
            V.append(V[-1] + new_V)
            S[-1] = S[-1] - new_V
            case = 2.1
        else:
            V.append(V[-1] + S[-1])
            S[-1] = 0
            case = 2.2

    elif t[-1] > (start + duration):
        V.append(V[-1])
        case = 3

    V_end = V[-1]
    if abs(V_start - V_end) > 50:
        print("Anomalous increase")
        print("Case: ", case)

meta_t = []
meta_S = []
meta_I = []
meta_R = []
meta_V = []
meta_T = []
runs = 20


for q in range(0,runs):
    N = 1000

    t = [0]
    S = [N]
    I = [1]
    R = [0]
    V = [0]
    T = [N]

    for i in range(0,8000):
        while t[-1] < 365*4:
            #Uniformly accelerate simulation
            acceleration = 4.9
            #Infection parameter
            beta = 0.030*acceleration
            #Recovery parameter
            gamma = 7*acceleration
            #Natural immunity loss parameter
            eta = 2*acceleration
            #Vaccination immunity loss parameter
            phi = 1.5*acceleration

            r1 = beta*S[-1]*I[-1]/N
            r2 = gamma*I[-1]/N
            r3 = eta*R[-1]/N
            r4 = eta*V[-1]/N

            big_R = r1 + r2 + r3 + r4

            if big_R != 0:
                dt = (1/big_R)*log(1/uni())

                #Generate the Monte Carlo random controller
                montepython = uni()*big_R

                #Susceptible person infected
                if montepython < r1:
                    t.append(t[-1] + dt)
                    S.append(S[-1] - 1)
                    I.append(I[-1] + 1)
                    R.append(R[-1])
                    V_population(S, V, t)
                    #Keep track of total population
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

                #Infected person recovered
                elif montepython < (r1 + r2):
                    t.append(t[-1] + dt)
                    S.append(S[-1])
                    I.append(I[-1] - 1)
                    R.append(R[-1] + 1)
                    V_population(S, V, t)
                    #Keep track of total population
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

                #Immune person loses immunity
                elif montepython < (r1 + r2 + r3):
                    t.append(t[-1] + dt)
                    S.append(S[-1] + 1)
                    I.append(I[-1])
                    R.append(R[-1] - 1)
                    V_population(S, V, t)
                    #Keep track of total population
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

                #Vaccinated person loses immunity
                elif V[-1] > 0:
                    t.append(t[-1] + dt)
                    S.append(S[-1] + 1)
                    I.append(I[-1])
                    R.append(R[-1])
                    V.append(V[-1] - 1)
                    #Keep track of total population
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

            else:
                break

    #Check data structure lengths for compatibility issues
    #print("t: ", len(t))
    #print("S: ", len(S))
    #print("I: ", len(I))
    #print("R: ", len(R))
    #print("V: ", len(V))

    #Save all realizations
    meta_t.append(t)
    meta_S.append(S)
    meta_I.append(I)
    meta_R.append(R)
    meta_V.append(V)
    meta_T.append(T)

    #Plot intermediate realization with some transparency
    plt.plot(t, S, color = 'green', alpha = 0.05)
    plt.plot(t, I, color = 'red', alpha = 0.05)
    plt.plot(t, list(map(add, R, V)), color = 'cyan', alpha = 0.05)
    plt.plot(t, R, '--', color = 'blue', alpha = 0.05)
    plt.plot(t, V, color = 'purple', alpha = 0.05)
    plt.plot(t, T, '--', color = 'black', alpha = 0.05)
    plt.title('SIRV Model')
    plt.xlabel("Time")
    plt.ylabel("People")

plt.plot(t, S, color = 'green', label = 'Susceptible')
plt.plot(t, I, color = 'red', label = 'Infected')
plt.plot(t, list(map(add, R, V))  , color = 'cyan', label = 'Immune')
plt.plot(t, R, '--', color = 'blue', label = 'Naturally Recovered')
plt.plot(t, V, color = 'purple', label = 'Vaccinated')
plt.plot(t, T, '--', color = 'black', label = 'Total')
plt.legend()
plt.show()

eradicated = 0
short = 0
short_criteria = 10

for i in range(0,runs):
    if len(meta_t[i]) < short_criteria:
        short += 1
    if len(meta_t[i]) > short_criteria and meta_I[i][-1] == 0:
        eradicated += 1

print("Statistics: ")
print("Pandemics that immediately fizzled: %f%%" % (100*short/runs))
print("Pandemics that resulted in eradication: %f%%" % (100*eradicated/(runs-short)))
print("Pandemics that continued: %f%%" % (100*( 1 - short/runs - eradicated/(runs-short) )))
