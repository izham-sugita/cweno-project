import numpy as np
import matplotlib.pyplot as plt

print(plt.rcParams.get('figure.figsize'))
pi = np.pi;
print(pi)

#Changing the default size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 20
fig_size[1] = 16
plt.rcParams["figure.figsize"] = fig_size

N = int(input("Enter maximum nodes\n"))

x1 = np.linspace(0.0, 2.0*pi,N)
x2 = np.linspace(0.0, 2.0*pi,N)

y1 = np.sin(x1)
y2 = np.sin(x1)
y3 = np.sin(x1)

imax = len(y1)
dx = x1[1]-x1[0]

print(dx)
print(len(y1))

for i in range(2,imax-3):
    x2[i] = x1[i] + 0.5*dx


for i in range(2,imax-3):
    y3[i] = 0.5*(y1[i+1] + y1[i])
    y2[i] = (-y1[i-2] + 7.0*y1[i] + 7.0*y1[i+1] -y1[i+2])/12.0
    #print(x1[i], x2[i])


'''    
plt.subplot(2, 1, 1)
plt.plot(x1, y1,'o-')
plt.plot(x2, y2, 'o')
plt.title('Interpolation comparison')
plt.ylabel('Y')

plt.subplot(2, 1, 2)
plt.plot(x2, y2, 'o')
plt.title('Interpolation data')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
'''
    
plt.plot(x1, y1, 'o-')
plt.plot(x2, y2, 'o')
plt.plot(x2, y3, '^')
plt.title('Interpolation sample')
plt.ylabel('y')
plt.show()

