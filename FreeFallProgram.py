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


y_init = 1000  # Initial Height m
vy_init = 0  # Initial velocity in ms^-1
g = 9.81  # Accel due to gravity ms^-2
m = 73  # Mass of object in kg
cd = 1.15 # Drag Coeff for a sky diver (Average)
rho0 = 1.2 # Air density at ambient temp + pressure in kg m^-3
h = 7640 # height at which rho = rho0 / e in m

# Number of calculations in each part
numOfPoints = 2000
# Number of dts in the time array
tStart = 0  # Starting time t
tEnd = 200 # End time t

# Taking Felix to be roughly the shape of an ellipse, where at largest area = 2pi*a*b
# His height, given by https://bodysize.org/en/felix-baumgartner/ is:
a = 1.7
# Estimated waist diameter, given by https://www.cdc.gov/nchs/fastats/mens-health.htm
b = 0.996 / np.pi
# So, his cross-sectional area is:
area = np.pi * a * b
# Constant drag factor, in kg m^-1
kConstant = cd * rho0 * area / 2


def partA(terminateAtFloor=True, variableRho=False):
    # Arrays holding time, displacement, and velocity in the y.
    global t
    global y
    global vy
    t = np.linspace(tStart, tEnd, numOfPoints)
    y = np.zeros(numOfPoints)
    vy = np.zeros(numOfPoints)

    # Applying initial value for y and vy:
    y[0] = y_init
    vy[0] = vy_init

    # Calculating using given formulae
    for i in range(1, numOfPoints):
        k = kConstant
        # If Rho is said to be variable, it will calculate it:
        if variableRho:
            k = (cd * area * rho0 / 2) * np.exp(-y[i-1]/h)

        vy[i] = -np.sqrt(m * g / k) * np.tanh(np.sqrt(k * g / m) * t[i - 1])
        y[i] = y_init - (m / k) * np.log(np.cosh(np.sqrt(k * g / m) * t[i - 1]))

        # If it goes below the floor, it will quit.
        if terminateAtFloor and y[i] <= 0:
            y = y[:i]
            vy = vy[:i]
            t = t[:i]
            break

    # Plotting
    plt.figure(1)
    plt.subplot(211)
    plt.title("Part A")
    plotData(t, y, 'Altitude /m', colour='red')
    plt.subplot(212)
    plotData(t, vy, 'Velocity /ms^-1', colour='pink')
    plt.show()


def partB(variableRho=False):
    # Arrays holding time, displacement, and velocity in the y.
    t = np.zeros(numOfPoints)
    y = np.zeros(numOfPoints)
    vy = np.zeros(numOfPoints)

    # Applying initial value for y and vy:
    y[0] = y_init
    vy[0] = vy_init

    dt = (tEnd-tStart)/numOfPoints

    for i in range(1, numOfPoints):
        k = kConstant
        # If Rho is said to be variable, it will calculate it:
        if variableRho:
            k = (cd * area * rho0 / 2) * np.exp(-y[i - 1] / h)

        vy[i] = vy[i-1] - dt*(g + (k / m) * np.abs(vy[i - 1]) * vy[i - 1])
        y[i] = y[i-1] + dt*(vy[i-1])
        t[i] = t[i-1] + dt

        # If it goes below the floor, it will quit.
        if y[i] <= 0:
            y = y[:i]
            vy = vy[:i]
            t = t[:i]
            break


    # Plotting
    plt.figure(2)
    plt.subplot(211)
    plt.title("Part B")
    plotData(t, y, 'bAltitude /m', colour='red')
    plt.subplot(212)
    plotData(t, vy, 'Velocity /ms^-1', colour='pink')
    plt.show()







# %%

print('------\nWelcome to Exercise 2\n\nObjectives:\n'
      '--Familiarity with Numpy and Scipy with ODES\n--Using Euler '
      'method for a 1D free fall problem including air resistance.')

while True:
    answer = input('\nPlease enter one of the following letters, or q to quit:\n --(a) : '
                   'Free fall with constant m, k, and g\n'
                   ' --(b) : Euler method for computing freefall\n'
                   ' --(c) : Euler method or Analytical method with variable air density\n')

    if answer == 'a':
        stopAtFloor = input("Do you want it to stop when it reaches the floor? Enter Y if so\n")
        partA(stopAtFloor == 'Y')

    elif answer == 'b':
        partB()

    elif answer == 'c':
        whichOne = input("Which of a or b would you like to see with variable air resistance?\n")
        if whichOne == 'a':
            partA(True, True)
        elif whichOne == 'b':
            partB(True)
        else:
            print("That was not an option, returning to start\n\n")

    else:
        print("Thank you")
        break
