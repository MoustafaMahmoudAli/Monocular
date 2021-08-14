import numpy as np
import matplotlib.pyplot as plt

#opern image
img = plt.imread("F:\\ENG\\MONO ROBO\\test\\samples\\horse.png")
ylim, xlim, dummy = img.shape
fig, ax = plt.subplots()
ax.imshow(img, extent=[0,xlim,ylim,0])

#read polygons
x = []
y = []
for t in open('polys.txt').read().split():
    a, b = t.strip('()').split(',')
    x.append(int(a))
    y.append(int(b))
x = np.array(x)
y = np.array(y)

#polt polygons over image
for i in range(0, len(x), 2):
    plt.plot(x[i:i+2], y[i:i+2],'r.-')

plt.xlim([-5,xlim+5])
plt.ylim([ylim+5,-5])
plt.show()
