import math
import pandas as pd
import matplotlib.pyplot as plt
#не забываем отредактировать working directory
data1 = pd.read_csv("plot1.csv")
data2 = pd.read_csv("plot2.csv")
data3 = pd.read_csv("plot3.csv")
data4 = pd.read_csv("plot4.csv")
X = [data1['ratio'], data2['ratio'], data3['ratio'], data4['ratio']] #первая колонка с абсциссами (отношение сигнал / шум)
Y = [data1['errors'], data2['errors'], data3['errors'], data4['errors']] #вторая колонка с ординатами (отношение ошиок к общему числу декодирований)
for i in range(len(X)):
    X[i] = [x for x in X[i] if Y[i][int(x * 10)] != 0]
    Y[i] = [y for y in Y[i] if y != 0]
    for j in range(len(Y[i])):
        Y[i][j] = math.log(Y[i][j], 10)

plt.plot(X[0], Y[0], 'r', label = 'size = 1') #строим график
plt.plot(X[3], Y[3], 'y', label = 'size = 2')
plt.plot(X[1], Y[1], 'b', label = 'size = 4')
plt.plot(X[2], Y[2], 'g', label = 'size = 16')

plt.legend()
plt.xlabel('Noise to Signal') #подписываем оси
plt.ylabel('Error to Total')
plt.title('Results')
plt.show() #показываем график