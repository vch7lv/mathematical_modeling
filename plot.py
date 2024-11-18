import matplotlib.pyplot as plt

allx = []
ally = []


with open('data.txt', 'r') as file:
    for line in file:
        x,y = map(float, line.split())
        allx.append(x)
        ally.append(y)

plt.plot(allx, ally)
plt.show()

