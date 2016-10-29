import sys
sys.path.append("/Library/Frameworks/Python.framework/Versions/3.5/lib/python3.5/site-packages")
import os
os.chdir("/Users/SagarJaiswal/Desktop/progammingShit")
import pylab as p
import math as m
massOfBall = 0.024
massOfDisc = 0.010
Diameter = 0.04
AirViscosity = 1.48e-5
springConstant = 490
AirDensity = 1.2
frontalArea = m.pi * ((Diameter/2)**2)
SpringFullLength = 100  # This is in millimeters, actual = 101.6
SpringMaxCompressionLength = 13  # This is in milimeqters, actual = 13.5
PVCCoefficientOfFriction = 0.4
NormalForce = ((massOfDisc+massOfBall)*9.81)
MaximumHeight = 2.0


def DetermineCombination(startHeight, Distance):
    angle = []
    changeInLength = []
    for z in range(SpringMaxCompressionLength, SpringFullLength):
        changeLen = SpringFullLength - z
        changeInLength.append((changeLen*(10**-3)))
    for i in range(0, 90):
        angleRad = i*(m.pi/180.0)
        angle.append(angleRad)
    labelNo = 0
    count = 0
    global Both
    Both = False
    while Both is False and count < len(angle):
        for y in changeInLength:
            DistanceCalc, SphereHeight = CalculateTraj(startHeight, y, Distance, angle[count])
            if DistanceCalc > Distance and DistanceCalc < Distance + 0.1 and SphereHeight > 0.44 and SphereHeight < 0.46:
                labelNo += 1
                Both = True
                angleFound = angle[count]
                angleFound1 = angleFound * (180.0/m.pi)
                changeInLengthSp = y
                compressionLength = (((changeInLengthSp*1000) - SpringFullLength)* -1)
                print("{}. Angle {} degrees, with the compressed length of {}mm".format(labelNo, angleFound1, compressionLength))
                CalculateTraj(startHeight, changeInLengthSp, Distance, angleFound, plot=True)
        count = count + 1


def ballDrag(Re):
    cd = (24/Re) + (2.6*(Re/5))/(1+(Re/5)**1.52)
    cd2 = cd + (0.411*(Re/263000)**-7.94)/(1+(Re/263000)**-8)
    return (cd2 + ((Re**0.8)/461000))


def InitialVelocity(springComp):
    TimeStep = 0.00001
    arrayLength = round(0.1/TimeStep)
    springLen = p.zeros(arrayLength)
    springLen[0] = springComp
    velocity = p.zeros(arrayLength)
    RE = p.zeros(arrayLength)
    ballDrag1 = p.zeros(arrayLength)
    totalForce = p.zeros(arrayLength)
    velocity[0] = 0.001
    RE[0] = (velocity[0] * Diameter)/AirViscosity
    ballDrag1[0] = ballDrag(RE[0])
    FrictionalForce = PVCCoefficientOfFriction * NormalForce * -1
    i = 0
    while (springLen[i] > 0):
        springForce = springLen[i] * springConstant
        springForce = springForce * 0.90
        DragForce = -0.5 * AirDensity * ((velocity[i])**2) * frontalArea * ballDrag1[i]
        totalForce[i] = springForce + DragForce + FrictionalForce
        acceleration = totalForce[i]/massOfBall
        velocity[i+1] = velocity[i] + (acceleration * TimeStep)
        springLen[i+1] = springLen[i] - (velocity[i] * TimeStep)
        RE[i] = velocity[i+1] * frontalArea/AirViscosity
        ballDrag1[i+1] = ballDrag(RE[i])
        i = i + 1
    return velocity[i]


def CalculateTraj(startHeight, springComp, Distance, Angle, plot=False):
    TimeStep = 0.001
    arrayLength = round(3/TimeStep)
    SphereHeight = p.zeros(arrayLength)
    SphereDistance = p.zeros(arrayLength)
    velocityX = p.zeros(arrayLength)
    velocityY = p.zeros(arrayLength)
    RE = p.zeros(arrayLength)
    ballDrag1 = p.zeros(arrayLength)
    SphereHeight[0] = startHeight
    SphereDistance[0] = 0
    inVelocity = InitialVelocity(springComp)
    velocityX[0] = inVelocity * m.cos(Angle)
    velocityY[0] = inVelocity * m.sin(Angle)
    RE[0] = m.sqrt((velocityX[0])**2 + (velocityY[0])**2) * (Diameter/AirViscosity)
    ballDrag1[0] = ballDrag(RE[0])
    target = 0
    i = 0
    angle = Angle * (180.0/m.pi)
    while (SphereDistance[i] < Distance+0.1 and SphereHeight[i] >= target and SphereHeight[i] < MaximumHeight and angle >= 25.0):
        velocityY[i+1] = velocityY[i] - (0.5*AirDensity*((velocityY[i])**2) * frontalArea * (ballDrag1[i]/massOfBall)) * TimeStep
        velocityX[i+1] = velocityX[i] - (0.5*AirDensity*((velocityX[i])**2) * frontalArea * (ballDrag1[i]/massOfBall)) * TimeStep
        velocityY[i+1] = velocityY[i+1] - (9.81*TimeStep)
        SphereDistance[i+1] = SphereDistance[i] + velocityX[i+1] * TimeStep
        SphereHeight[i+1] = SphereHeight[i] + velocityY[i+1] * TimeStep
        if SphereHeight[i] > startHeight:
            if SphereHeight[i] > 0.45:
                target = 0.45
            else:
                target = SphereHeight[i]
        RE[i+1] = (m.sqrt((velocityX[i+1])**2 + (velocityY[i+1])**2) * (Diameter/AirViscosity))
        ballDrag1[i+1] = ballDrag(RE[i+1])
        i = i + 1
    if SphereHeight[i] < 0.45:
        global Both
        Both = False
    if plot is True:
        plot = False
        plotGraph(SphereDistance, SphereHeight, i)
    return SphereDistance[i], SphereHeight[i]


def plotGraph(x, y, i):
    p.plot(x[0:i], y[0:i])
