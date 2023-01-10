from random import uniform
from operator import add
from math import log
import matplotlib.pyplot as plt


def uni():
    #Input: None.
    #Output: Random float from 0 to 1.
    #Debugging values: None.

    #Allows user to map to map the uniform probability distribution to any
    #interval from 0 to t by mutliplying the uni() instance of this method
    #by the parameter t.
    return uniform(0,1)

def V_population(S, V, t):
    #Input: Record of the susceptible and vaccinated populations, and the
    #time intervals.
    #Output: Possible output of debugging values.
    #Debugging values: Case variable. If vaccination numbers change rapidly,
    #the regime of the vaccination period is recorded and outputed.
    #There have been previous bugs with the vaccinated population.

    #Debugging variable.
    V_start = V[-1]
    #Debugging variable.
    case = 0

    #Percentage of final vaccinated people.
    vax = 0.3
    #Start date of vaccination campaing.
    start = 360
    #Duration of vaccination campaign.
    duration = 365

    #If the last time is still before the vaccination campaign period.
    if t[-1] < start:
        V.append(0)
        #Debugging variable.
        case = 1

    #If the last time is during the vaccination campaign period.
    elif t[-1] >= start and t[-1] < (start + duration):
        #How many days have passed since the beginning of the campaign.
        t0 = t[-1] - start
        #How many people should be vaccinated by now. Equivalent to making
        #use of all vaccine doses used so far.
        expected_V = N*vax*t0/duration
        #How many unused doses are there. Produced vaccines minus vaccinated people.
        #Eventually replace V[-1] with a count of vaccines applied so far because
        #the V array is losing population to the progressive decay in immunity.
        new_V = int(expected_V-V[-1])

        #If there are enough susceptible people to apply all the new_V doses.
        if S[-1] > new_V:
            #Add new_V to the V population and subtract it from the
            #susceptible population.
            V.append(V[-1] + new_V)
            #The susceptible population is edited because it was saved before
            #calling the V_population method.
            S[-1] = S[-1] - new_V
            #Debugging variable.
            case = 2.1
        else:
            #If there are more vaccines available than susceptible population,
            #the susceptible population is added to the vaccinated population
            #and the former becomes 0.
            V.append(V[-1] + S[-1])
            S[-1] = 0
            #Debugging variable.
            case = 2.2

    #If the last time is after the vaccination campaign period, vaccinated
    #population remains the same.
    elif t[-1] > (start + duration):
        V.append(V[-1])
        #Debugging variable.
        case = 3

    V_end = V[-1]
    #Check for anomalous changes in vaccinated population numbers within this
    #instance of the V_population method.
    if abs(V_start - V_end) > 50:
        print("Anomalous change in V")
        print("Case: ", case)

#Memory of realization state variables and points in time.
meta_t = []
meta_S = []
meta_I = []
meta_R = []
meta_V = []
meta_T = []

#Number of realizations.
runs = 20


for q in range(0,runs):
    #Total initial population.
    N = 1000

    t = [0]
    #Disease begins with one person.
    I = [1]
    #SUsceptiuble people are the whole population minus the initial infected population.
    S = [N - I[-1]]
    R = [0]
    V = [0]
    #Total population tracker for debugging.
    T = [N]

    #Generate, at most 8000 time points per run.
    for i in range(0,8000):
        #Limit time to 4 'years' to introduce a standard time to simulations
        #and avoid having to compare 4 day realizations to 20 year realizations.
        while t[-1] < 365*4:
            #Uniformly increase all parameters and 'scale' the model in time so
            #cyclical events and the final, stable disease can be observed in
            #the four year window. One should not compare models with different
            #acceleration parameters.
            acceleration = 4.9
            #Infection parameter.
            beta = 0.030*acceleration
            #Recovery parameter.
            gamma = 7*acceleration
            #Natural immunity loss parameter.
            eta = 2*acceleration
            #Vaccination immunity loss parameter.
            phi = 1.5*acceleration

            #New infections are proportional to S and I populations.
            r1 = beta*S[-1]*I[-1]/N
            #Recuperations are proportional to I population.
            r2 = gamma*I[-1]/N
            #Natural immunity loss is proportional to R population.
            r3 = eta*R[-1]/N
            #Vaccination immunity loss is proportional to V population.
            r4 = phi*V[-1]/N

            big_R = r1 + r2 + r3 + r4

            if big_R != 0:
                #Generate poisson distributed time interval based on the
                #mean likelihood of the ocurrence of an event.
                dt = (1/big_R)*log(1/uni())

                #Generate the Monte Carlo random variable.
                montepython = uni()*big_R

                #Susceptible person infected.
                if montepython < r1:
                    t.append(t[-1] + dt)
                    S.append(S[-1] - 1)
                    I.append(I[-1] + 1)
                    R.append(R[-1])
                    #Check vaccinated population.
                    V_population(S, V, t)
                    #Keep track of total population.
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

                #Infected person recovered.
                elif montepython < (r1 + r2):
                    t.append(t[-1] + dt)
                    S.append(S[-1])
                    I.append(I[-1] - 1)
                    R.append(R[-1] + 1)
                    #Check vaccinated population.
                    V_population(S, V, t)
                    #Keep track of total population.
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

                #Immune person loses immunity.
                elif montepython < (r1 + r2 + r3):
                    t.append(t[-1] + dt)
                    S.append(S[-1] + 1)
                    I.append(I[-1])
                    R.append(R[-1] - 1)
                    #Check vaccinated population.
                    V_population(S, V, t)
                    #Keep track of total population.
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

                #Vaccinated person loses immunity.
                elif V[-1] > 0:
                    t.append(t[-1] + dt)
                    S.append(S[-1] + 1)
                    I.append(I[-1])
                    R.append(R[-1])
                    #Vaccinated population is not changed because it would
                    #revert the loss just computed.
                    V.append(V[-1] - 1)
                    #Keep track of total population.
                    T.append(S[-1] + I[-1] + R[-1] + V[-1])

            #If time exceeds 'four years', simulation ends.
            else:
                break

    #Check data structure lengths for compatibility issues.
    #print("t: ", len(t))
    #print("S: ", len(S))
    #print("I: ", len(I))
    #print("R: ", len(R))
    #print("V: ", len(V))

    #Save realization to memory.
    meta_t.append(t)
    meta_S.append(S)
    meta_I.append(I)
    meta_R.append(R)
    meta_V.append(V)
    meta_T.append(T)

    #Plot the realization with transparency to plt.
    plt.plot(t, S, color = 'green', alpha = 0.05)
    plt.plot(t, I, color = 'red', alpha = 0.05)
    plt.plot(t, list(map(add, R, V)), color = 'cyan', alpha = 0.05)
    plt.plot(t, R, '--', color = 'blue', alpha = 0.05)
    plt.plot(t, V, color = 'purple', alpha = 0.05)
    plt.plot(t, T, '--', color = 'black', alpha = 0.05)
    plt.title('SIRV Model')
    plt.xlabel("Time")
    plt.ylabel("People")

#Plot the last realization to the plt. Eventually must change this to an average
#of the realizations. The step after that is to plot the average of the endemic
#disease realizations and eradication realizations separetely.
plt.plot(t, S, color = 'green', label = 'Susceptible')
plt.plot(t, I, color = 'red', label = 'Infected')
plt.plot(t, list(map(add, R, V))  , color = 'cyan', label = 'Immune')
plt.plot(t, R, '--', color = 'blue', label = 'Naturally Recovered')
plt.plot(t, V, color = 'purple', label = 'Vaccinated')
plt.plot(t, T, '--', color = 'black', label = 'Total')
plt.legend()
plt.show()

#Initialize the variables for counting eradication scenarios and endemic
#disease scenarios, and the maximum days for a disease to be considered to have
#'fizzledout'.
short_criteria = 10
short = 0
eradicated = 0

for i in range(0,runs):
    if len(meta_t[i]) < short_criteria:
        short += 1
    if len(meta_t[i]) > short_criteria and meta_I[i][-1] == 0:
        eradicated += 1

print("Statistics: ")
#Short scenarios out of all runs.
print("Pandemics that immediately fizzled: %f%%" % (100*short/runs))
#Eradication scenarios out of all runs that did not fizzle out.
print("Pandemics that resulted in eradication: %f%%" % (100* eradicated/(runs-short) ))
#Endemic disease scenarios out of all runs that did not fizzle out.
print("Pandemics that continued: %f%%" % (100* (runs-short-eradicated)/(runs-short) ))
