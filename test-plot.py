import numpy as np
import matplotlib.pyplot as plt

def test(number):
    for i in range(len(number)):
        print(number[i])
        

def cweno4(stencil):
    um2 = stencil[0]
    um1 = stencil[1]
    u0 = stencil[2]
    up1 = stencil[3]
    up2 = stencil[4]

    eps = 1.0e-6
    c0 = 0.16666666 #1.0/6.0
    c1 = 0.66666666 #2.0/3.0
    c2 = 0.16666666 #1.0/6.0

    #13/12 = 1.08333333

    b0 = 0.25*(um2-4.0*um1+3.0*u0)**2 + 1.08333333*(um2-2.0*um1+u0)**2
    b1 = 0.25*(up1-um1)**2 + 1.08333333*(um1-2.0*u0+up1)**2
    b2 = 0.25*(3.0*u0-4.0*up1+up2)**2 + 1.08333333*(u0-2.0*up1+up2)**2

    a0 = c0/( ( eps + b0)**2  )
    a1 = c1/( ( eps + b1)**2  )
    a2 = c2/( ( eps + b2)**2  )

    w0 = a0/( a0+a1+a2  )
    w1 = a1/( a0+a1+a2  )
    w2 = a2/( a0+a1+a2  )

    uhalf = 0.5*( -w0*um1 + (3.0*w0+w1)*u0 +(3.0*w2+w1)*up1 -w2*up2  )
    
    return uhalf


def cweno4_int(stencil):
    interpolate = 1.0

    return interpolate

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
stencil = np.ndarray((5), dtype=float)


y1 = np.sin(x1)
y2 = np.sin(x1)
y3 = np.sin(x1)
y4 = np.sin(x1)

'''
for i in range(N):
    y1[i] = 0.0
    y4[i] = 0.0
    if(x1[i] > 0.5*pi and x1[i] < pi):
        y1[i] = 1.0
'''        

f = open("init.csv", "w")
f.write("x, y\n")
for i in range(N):
    xg = str(x1[i])
    yg = str(y1[i])
    st0 = xg+","+yg+"\n"
    f.write(st0)
    
f.close()

'''
plt.plot(x1, y1, 'ko') #analytical 
plt.title('Function initial condition')
plt.ylabel('y')
plt.show()
'''

'''
y2 = y1
y3 = y1
y4 = y1
'''

imax = len(y1)
dx = x1[1]-x1[0]

print(dx)
print(len(y1))

for i in range(2,imax-3):
    y3[i] = 0.5*(y1[i+1] + y1[i])
    y2[i] = (-y1[i-2] + 7.0*y1[i]
             + 7.0*y1[i+1] -y1[i+2])/12.0 #central interpolation
    err1 = abs(y3[i]-y1[i])
    err2 = abs(y2[i]-y1[i])
    #print("Ave. based error=%f, 4th order based=%f"%(err1, err2))


for i in range(3, imax-4):
    '''
    #for mid-point derivative reconstruction
    stencil[0] = (y1[i-1] - y1[i-2])/dx
    stencil[1] = (y1[i] - y1[i-1])/dx
    stencil[2] = (y1[i+1] - y1[i])/dx
    stencil[3] = (y1[i+2] - y1[i+1])/dx
    stencil[4] = (y1[i+3] - y1[i+2])/dx 
    '''  

    #for mid-point value reconstruction
    stencil[0] = y1[i-2]
    stencil[1] = y1[i-1]
    stencil[2] = y1[i]
    stencil[3] = y1[i+1]
    stencil[4] = y1[i+2]
    

    #interpolation
    y4[i] = cweno4(stencil)
    x2[i] = x1[i] + 0.5*dx

    #error = abs(y4[i] - np.cos(x2[i]))
    error = abs( y4[i]-np.sin(x2[i]) )
    print(error)


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

'''
import matplotlib.lines as mlines
blue_line = mlines.Line2D([], [], color='blue', marker='o',
                          markersize=15, label='Analytical')
plt.legend(handles=[blue_line])
'''

'''
plt.plot(x1, y1, 'ko-') #analytical 
plt.plot(x2, y4, 'ro') #4th order central WENO interpolation
#plt.plot(x2, y2, 'b^') #4th order CD interpolation
#plt.plot(x2, y3, '^') #average [i,i+1]
plt.title('Interpolation sample')
plt.ylabel('y')
plt.show()
'''
