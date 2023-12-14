import re
import matplotlib.pyplot as plt
import numpy
import random
import sys
def f(x):
    return float(1 - numpy.exp(-x)) / (1 + numpy.exp(-x))

def fDerivative(x):
    return float((1 + f(x))*(1 - f(x))) / 2

def readFiles(address):
    fileName = address
    parameters = []
    with open(fileName,"r") as file :
        lines = file.readlines()
        temp = []
        help = []
        for line in lines :
            temp.append(re.split(r'\t|\n', line))
        for help in temp:
            help.remove("")
        for x in temp:
            parameters.append([float(i) for i in x]) 
    return parameters

def MLP(parameters, alpha, cellsNo):
    count = 0
    w = [None] * (cellsNo + 1)
    v = [[[None] for i in range(3)] for j in range(cellsNo)]
    y, y_in = 0, 0
    z, z_in = [None] * cellsNo, [None] * cellsNo
    deltaK, deltaW = 0, [None] * (cellsNo + 1)
    delta_inJ, deltaJ = [None] * cellsNo, [None] * cellsNo
    deltaV = [[[None] for i in range(3)] for j in range(cellsNo)]
    epochs, oldError, newError = 0, 0, 999999
    # random.shuffle(parameters)
    index = int(len(parameters) * 0.7)
    samples = parameters[:index]
    validations = parameters[index:]
    for i in range(cellsNo + 1):
        w[i] = random.uniform(-0.5, 0.5)
    for i in range(cellsNo):
        for j in range(3):
            v[i][j] = random.uniform(-0.5, 0.5)
    while(True):
        epochs += 1
        for point in samples:
            for i in range(cellsNo):
                z_in[i] = 0
                for j in range(1,3):
                    z_in[i] += point[j-1] * v[i][j]
                z_in[i] += v[i][0]
            for i in range(cellsNo):
                z[i] = f(z_in[i])
            y_in = 0
            for i in range(cellsNo):
                y_in += z[i] * w[i + 1]
            y_in += w[0]
            y = f(y_in)
            # ----------------------------------------------------------------
            deltaK = (point[2] - y) * fDerivative(y_in)
            for i in range(cellsNo):
                deltaW[i + 1] = alpha * deltaK * z[i]
            deltaW[0] = alpha * deltaK

            for i in range(cellsNo):
                delta_inJ[i] = deltaK * w[i + 1]
            for i in range(cellsNo):
                deltaJ[i] = delta_inJ[i] * fDerivative(z_in[i])

            for i in range(cellsNo):
                for j in range(1,3):
                    deltaV[i][j] = alpha * deltaJ[i] * point[j - 1]
                deltaV[i][0] = alpha * deltaJ[i]
            # ----------------------------------------------------------------
            for i in range(cellsNo + 1):
                w[i] += deltaW[i]
            for i in range(cellsNo):
                for j in range(3):
                    v[i][j] += deltaV[i][j]
        oldError = newError
        newError = 0
        for point in validations:
            for i in range(cellsNo):
                z_in[i] = 0
                for j in range(1,3):
                    z_in[i] += point[j-1] * v[i][j]
                z_in[i] += v[i][0]
            for i in range(cellsNo):
                z[i] = f(z_in[i])
            y_in = 0
            for i in range(cellsNo):
                y_in += z[i] * w[i + 1]
            y_in += w[0]
            y = f(y_in)
            e = float(point[2] - y)
            newError += (((e ** 2) / 2))
        if newError > oldError:
            count += 1
        else :
            count = 0
        if count == 10 :
            break
    weights = []
    weights.append(w)
    weights.append(v)
    weights.append(epochs)
    return weights

def position(cellsNo,w,v,point):
    z, z_in = [None] * cellsNo, [None] * cellsNo
    for i in range(cellsNo):
        z_in[i] = 0
        for j in range(1,3):
            z_in[i] += point[j - 1] * v[i][j]
        z_in[i] += v[i][0]
    for i in range(cellsNo):
        z[i] = f(z_in[i])
    y_in = 0
    for i in range(cellsNo):
        y_in += z[i] * w[i + 1]
    y_in += w[0]
    y = f(y_in)
    if y > 0 :
        return 2
    else:
        return -2

def main1():
    # for i in range(1,7):
        # parameters = readFiles(str(i) + ".txt")
        parameters = readFiles("5.txt")
        alpha = 0.0001
        cellsNo = 3
        weights = MLP(parameters, alpha, cellsNo)
        # print(weights[0])
        # print(weights[1])
        # print(weights[2])
        xlimits = [[-0.8,1.2] ,[-1,1.2] ,[-1,1] ,[-1,1.2] ,[-0.6,1] ,[-0.8,1.2]]
        ylimits = [[-1,1.2] ,[-1,1.2] ,[-0.8,1.2] ,[-1.2,0.6] ,[-1,1.2] ,[-1.2,1]]
        x1 ,x2 ,y1 ,y2 =([] for i in range(4))
        x3 , x4, y3, y4 = ([] for i in range(4))
        # for j in numpy.arange(xlimits[i-1][0],xlimits[i-1][1],0.1):
        #     for k in numpy.arange(ylimits[i-1][0],ylimits[i-1][0],0.1):
        for j in numpy.arange(xlimits[4][0],xlimits[4][1],0.05):
            for k in numpy.arange(ylimits[4][0],ylimits[4][1],0.05):
                answer = position(cellsNo, weights[0], weights[1], [j,k])
                if answer == 2:
                    x3.append(j)
                    y3.append(k)
                else:
                    x4.append(j)
                    y4.append(k)
        for x in parameters:
            if(x[2] == 1):
                x1.append(x[0])
                y1.append(x[1])
            else:
                x2.append(x[0])
                y2.append(x[1])

        # if(i == 3): X = [-0.5,0.5]
        # else: X = [-1,1]
        # X = [-1,1]
        # for j in range(cellsNo):
        #     plt.plot(X,[(-weights[1][j][1]*help -weights[1][j][0]+1)/weights[1][j][2] for help in X], linestyle = ':',c='#4CAF50')
        # plt.plot(X,[(-weights[0]*help -weights[2]+1)/weights[1] for help in X], linestyle = ':',c='#4CAF50')
        # plt.plot(X,[(-weights[0]*help -weights[2])/weights[1] for help in X])
        # plt.plot(X,[(-weights[0]*help -weights[2]-1)/weights[1] for help in X], linestyle = ':',c='#4CAF50')
        plt.xlabel("X1")
        plt.ylabel("X2")
        # plt.title("ّFile " + str(i) + "         Accuracy = " + str(accuracy))
        # plt.title("1.txt")
        # plt.grid(linestyle = '--', linewidth = 0.5)
        plt.scatter(x1,y1)
        plt.scatter(x2,y2)
        plt.scatter(x3,y3,s = 5)
        plt.scatter(x4,y4 , s = 5)
        # plt.xlim(xlimits[i-1])
        # plt.ylim(ylimits[i-1])
        plt.xlim(xlimits[4])
        plt.ylim(ylimits[4])
        plt.show() 
        

main1()