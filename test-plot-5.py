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
    c0 = 1.0/6.0
    c1 = 2.0/3.0
    c2 = 1.0/6.0

    #13/12 = 1.08333333

    b0 = 0.25*(um2-4.0*um1+3.0*u0)**2 + (13.0/12.0)*(um2-2.0*um1+u0)**2
    b1 = 0.25*(up1-um1)**2 + (13.0/12.0)*(um1-2.0*u0+up1)**2
    b2 = 0.25*(3.0*u0-4.0*up1+up2)**2 + (13.0/12.0)*(u0-2.0*up1+up2)**2

    a0 = c0/( ( eps + b0)  )
    a1 = c1/( ( eps + b1)  )
    a2 = c2/( ( eps + b2)  )

    w0 = a0/( a0+a1+a2  )
    w1 = a1/( a0+a1+a2  )
    w2 = a2/( a0+a1+a2  )

    uhalf = 0.5*( -w0*um1 + (3.0*w0+w1)*u0 +(3.0*w2+w1)*up1 -w2*up2  )
    
    return uhalf


def cweno4_my(stencil):
    um2 = stencil[0]
    um1 = stencil[1]
    u0 = stencil[2]
    up1 = stencil[3]
    up2 = stencil[4]

    eps = 1.0e-6
    c0 = 1.0/6.0
    c1 = 2.0/3.0
    c2 = 1.0/6.0

    p0 = 0.5*(3.0*u0 - 4.0*um1 + um2) # /dx comes from stencil input
    p1 = 0.5*(up1 - um1)
    p2 = 0.5*(-3.0*u0 + 4.0*up1 - up2)

    #13/12 = 1.08333333

    b0 = 0.25*(um2-4.0*um1+3.0*u0)**2 + (13.0/12.0)*(um2-2.0*um1+u0)**2
    b1 = 0.25*(up1-um1)**2 + (13.0/12.0)*(um1-2.0*u0+up1)**2
    b2 = 0.25*(3.0*u0-4.0*up1+up2)**2 + (13.0/12.0)*(u0-2.0*up1+up2)**2

    a0 = c0/( ( eps + b0)  )
    a1 = c1/( ( eps + b1)  )
    a2 = c2/( ( eps + b2)  )

    w0 = a0/( a0+a1+a2  )
    w1 = a1/( a0+a1+a2  )
    w2 = a2/( a0+a1+a2  )

    interpolate = w0*p0 + w1*p1 + w2*p2
    
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

dy4dx = np.cos(x1)

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

imax = len(y1)
dx = x1[1]-x1[0]

print(dx)
print(len(y1))


for i in range(3,imax-4):
    y2[i] = 0.5*(y1[i+1] + y1[i])
    
    #interpolation x_(i+1/2)
    stencil[0] = y1[i-2]
    stencil[1] = y1[i-1]
    stencil[2] = y1[i]  #interpolation for point [i,i+1]
    stencil[3] = y1[i+1]
    stencil[4] = y1[i+2]

    yplushalf = cweno4(stencil)
    
    y4[i] = yplushalf
    x2[i] = x1[i] + 0.5*dx

    
    #err1 = abs( y3[i] -np.cos(x1[i]) )
    #err2 = abs( y2[i] -np.cos(x1[i]) )
    #err3 = abs( y4[i] -np.cos(x2[i]) )


sum = 0.0
for i in range(3, imax-4):
#stencil for reconstruction

    #stencil[0] = y4[i-2]/dx
    #stencil[1] = y4[i-1]/dx
    #stencil[2] = y4[i]/dx
    #stencil[3] = y4[i+1]/dx
    #stencil[4] = y4[i+2]/dx

    #dy4dx[i] = cweno4_my(stencil)
    #err = (np.cos(x2[i]) - dy4dx[i])**2

    stencil[0] = (y4[i-2]-y4[i-3]) /dx
    stencil[1] = (y4[i-1]-y4[i-2]) /dx
    stencil[2] = (y4[i]-y4[i-1]) /dx #reconstruction [i-1,i,i+1] from the interpolation point
    stencil[3] = (y4[i+1] - y4[i]) /dx
    stencil[4] = (y4[i+2] - y4[i+1])/dx

    dy4dx[i] = cweno4(stencil)


f = open("derivatives.csv","w")
f.write("x, dy4dx, analytic\n")
for i in range(3, imax-4):
    str0 = str(x2[i])
    str1 = str(dy4dx[i])
    str2 = str(np.cos(x2[i]))
    strall = str0+","+str1+","+str2+"\n"
    f.write(strall)

f.close()

plt.plot(x2, y4, 'ko-') #analytical 
plt.plot(x2, dy4dx, 'r^') #4th order central WENO interpolation
plt.title('Interpolation sample')
plt.ylabel('y')
plt.show()


f =open("interpolation.csv","w")
f.write("x, 4thCWENO, 2ndCD, analytic\n")
for i in range(3, imax-4):
    str0 = str(x2[i])
    str1 = str(y4[i])
    str2 = str(y2[i])
    str3 = str(np.sin(x2[i]))
    strall = str0+","+str1+","+str2+","+str3+"\n"
    f.write(strall)
    #print(y4[i] - np.cos(x2[i]))

f.close()

    
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
