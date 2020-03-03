import numpy as np
import matplotlib.pyplot as plt

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

def cweno4_v2(stencil):

    um1 = stencil[0]
    u0 = stencil[1]
    up1 = stencil[2]
    up2 = stencil[3]

    eps = 1.0e-6
    c0 = 0.125 #1.0/8.0
    c1 = 0.75 #3.0/4.0
    c2 = 0.125 #1.0/8.0

    #13/12 = 1.08333333
    
    b0 = (u0 - um1)**2
    b1 = (up1 - u0)**2 
    b2 = (up2 - up1)**2 

    a0 = c0/( ( eps + b0)**2  )
    a1 = c1/( ( eps + b1)**2  )
    a2 = c2/( ( eps + b2)**2  )

    w0 = a0/( a0+a1+a2  )
    w1 = a1/( a0+a1+a2  )
    w2 = a2/( a0+a1+a2  )
    
    #w0 = c0
    #w1 = c1
    #w2 = c2

    uhalf = w0*(1.5*u0 - 0.5*um1 ) + 0.5*w1*(u0+up1) + w2*(1.5*up1 - 0.5*up2)  

    return uhalf


print(plt.rcParams.get('figure.figsize'))

#Changing the default size
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 20
fig_size[1] = 16
plt.rcParams["figure.figsize"] = fig_size

pi = np.pi;

#N = int(input("Enter maximum nodes\n"))
N = int(input("Enter the first iteration nodes\n"))

meshlist = [N, 2*(N-1)+1, 3*(N-1)+1, 4*(N-1)+1]

jmax = len(meshlist)
logdx = np.zeros(jmax, dtype=float)
logerr = np.zeros(jmax, dtype=float)

j = 0

for imax in meshlist:
    x1 = np.linspace(0.0, 2.0*pi,imax)
    x2 = np.linspace(0.0, 2.0*pi,imax)
    stencil = np.ndarray((5), dtype=float)
    dx = x1[1]-x1[0]

    y1 = np.sin(x1)
    y2 = np.sin(x1)
    y3 = np.cos(x1)
    y4 = np.cos(x1)

    for i in range(3, imax-4):
                
        #Interpolate mid-values
        stencil[0] = y1[i-2] 
        stencil[1] = y1[i-1] 
        stencil[2] = y1[i] 
        stencil[3] = y1[i+1] 
        stencil[4] = y1[i+2] 
        
        y2[i] = cweno4(stencil)
        #x2[i] = x1[i] + 0.5*dx
        #err = ( y2[i] - np.cos(x2[i]) )**2
        

    for i in range(3, imax-4):
        #Use interpolated mid-values to reconstruct derivatives at mid-point
        stencil[0] = (y2[i-2] -y2[i-3])/dx
        stencil[1] = (y2[i-1] - y2[i-2])/dx
        stencil[2] = (y2[i] - y2[i-1])/dx 
        stencil[3] = (y2[i+1] - y2[i])/dx 
        stencil[4] = (y2[i+2] -y2[i+1])/dx 
        
        y3[i] = cweno4(stencil)
        x2[i] = x1[i] + 0.5*dx
        err = ( y3[i] - np.cos(x2[i]) )**2

    err = np.sqrt(err/(N-7 ) )
    logdx[j] = abs(np.log(dx))
    logerr[j] = abs(np.log(err))
    j +=1
    
    print(imax,dx,err)
            

A = np.vstack( [logdx, np.ones(len(logdx)) ]  ).T
m, c = np.linalg.lstsq(A, logerr, rcond=None)[0]
print(m,c)

linelabel = str(m)+"x +"+str(c)
print(linelabel)

plt.plot(logdx, logerr, 'o', label='Original data', markersize=10)
plt.plot(logdx, m*logdx + c, 'r', label=linelabel)
plt.gca().invert_xaxis()
plt.legend()
plt.show()

'''
plt.plot(x1, y1, 'ko-') #analytical
plt.plot(x2, y2, 'b^') #interpolation 
plt.title('Function interpolation')
plt.ylabel('y')
plt.show()
'''
