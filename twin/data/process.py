#!/opt/local/bin/python3.8
import numpy as np
import sys
from math import sin,cos,pi,sqrt,atan2
from struct import pack

# Number of beads
n=2

# Starting frame
s=0

# Load tracking data and measure its size
f=np.load("2020_09_15/n_beads_2/n_beads_2_velocity_150_trajectories.npz")
l=len(f['Time'])//n

# Extract times
a=np.empty((l),dtype=np.float32)
b=f['Time']
c=f['ID']

j=0
while j<l:

    # Check all of the IDs are in order, and that the times agree
    if c[j*n]!=0:
        break
    stop=False
    for i in range(1,n):
        if c[j*n+i]!=i or b[j*n+i]!=b[j*n]:
            stop=True
            break

    # A full set of n beads has been found. Record the time.
    a[j]=b[j*n]
    j+=1

# Print status message about what fraction of the file has been read
print(j,"records out of a potential",l,"processed")
if j<=s:
    print("No records to save")
    sys.exit()

# Output the times
b=open("2_tfix.dat","wb")
b.write(pack('ii',n,j-s))
b.write(a[s:j].tobytes())

# Convert positions to single-precision floats
x=f['X'][s*n:j*n]+0.5
y=f['Y'][s*n:j*n]+0.5

for i in range((j-s)*n):
    print(x[i],y[i])

cx=0.5*(np.amax(x)+np.amin(x))
cy=0.5*(np.amax(y)+np.amin(y))

x[:]-=cx
y[:]-=cy
xs=0
ys=0
for i in range((j-s)*n):
    r=sqrt(x[i]*x[i]+y[i]*y[i])
    th=atan2(y[i],x[i])
    xs+=r*cos(6*th)
    ys+=r*sin(6*th)

ro=atan2(ys,xs)/6.
print(ro*180/pi)

cr=cos(ro)
sr=sin(ro)
for i in range((j-s)*n):
    xx=cx+cr*x[i]-sr*y[i]
    yy=cy+sr*x[i]+cr*y[i]
    x[i]=xx
    y[i]=yy

print(cx,cy)
print(np.amin(x),np.amax(x))
print(np.amin(y),np.amax(y))

xy=np.empty(((j-s)*n,2),dtype=np.float32)
xy[:,0]=x
xy[:,1]=y

for i in range(n):
    print(xy[i,0],xy[i,1])

print("\n\n")
for i in range(n):
    print(xy[i+n,0],xy[i+n,1])

print("\n\n")
for i in range(n):
    print(xy[i+2*n,0],xy[i+2*n,1])

b.write(xy.tobytes())

c=open("2020_09_15/n_beads_2/n_beads_2_velocity_150_theta_time.txt","r")

for i in range(s+1):
    c.readline()

z=np.empty((j-s),dtype=np.float32)
for i in range(j-s):
    a=c.readline().split()
    if len(a)<2:
        print("Terminated with",i,"of",j-s)
        break
    z[i]=float(a[1])
b.write(z.tobytes())

b.close()
