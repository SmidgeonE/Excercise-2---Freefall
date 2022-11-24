# %%

# Imports
import numpy as np
import matplotlib.pyplot as plt


# %%

# Reused function for plotting, for brevity
def plotData(xData, yData, yLabel, xLabel='Time /s', colour='green'):
    plt.plot(xData, yData, color=colour)
    plt.ylabel(yLabel)
    plt.xlabel(xLabel)


def partA(terminateAtFloor=True):
    y_init = 5000  # Initial Height m
    vy_init = 0  # Initial velocity in ms^-1
    g = 9.81  # Accel due to gravity ms^-2
    m = 73  # Mass of object in kg
    k = 10  # Constant drag factor, in kg m^-1

    numOfPoints = 200  # Number of dts in the time array
    tStart = 0  # Starting time t
    tEnd = 10  # End time t

    # Arrays holding time, displacement, and velocity in the y.
    t = np.linspace(tStart, tEnd, numOfPoints)
    y = np.zeros(numOfPoints)
    vy = np.zeros(numOfPoints)

    # Applying initial value for y and vy:
    y[0] = y_init
    vy[0] = vy_init

    # Calculating using given formulae
    for i in range(1, numOfPoints):
        vy[i] = -np.sqrt(m * g / k) * np.tanh(np.sqrt(k * g / m) * t[i - 1])
        y[i] = y[i - 1] - (m / k) * np.log(np.cosh(np.sqrt(k * g / m) * t[i - 1]))

        # If it goes below the floor, it will quit.
        if terminateAtFloor and y[i] <= 0:
            y = y[:i]
            vy = vy[:i]
            break

    # Plotting
    plt.figure(1)
    plt.subplot(211)
    plotData(t, y, 'Altitude /m', colour='red')
    plt.subplot(212)
    plotData(t, vy, 'Velocity /ms^-1', colour='pink')
    plt.show()


# %%

print('------\nWelcome to Exercise 2\n\nObjectives:\n'
      '--Familiarity with Numpy and Scipy with ODES\n--Using Euler '
      'method for a 1D free fall problem including air resistance')

answer = input('\nPlease enter one of the following letters, or q to quit:\n --(a) : '
               'Free fall with constant m, k, and \n')
if answer == 'a':
    partA(False)
