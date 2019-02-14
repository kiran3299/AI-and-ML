import numpy as np 
import matplotlib.pyplot as plt
import math 

A=np.array([[0,1],[0,-3]])
B=np.array([[5,0],[0,1]])
C=np.array([[0,-1],[2,0]])

P=B-A
Q=C-A

T=P[0,:]
Y=Q[1,:]
T1=P[1,:]
Y1=Q[0,:]

Temp1=0
Temp2=0
Temp3=0

for i in range(0,2):
	for j in range(0,2):
		if (i+j==0):
			Temp1=Temp1+T[i]*Y[j]-T1[i]*Y1[j]
		elif (i+j==1):
			Temp2=Temp2+T[i]*Y[j]-T1[i]*Y1[j]
		elif (i+j==2):
			Temp3=Temp3+T[i]*Y[j]-T1[i]*Y1[j]
			
	
a=Temp3
b=Temp2
c=Temp1-56

k1=(-b+ (math.sqrt(((b)**2)-4.0*a*(c))))/(2.0*a)
k2=(-b- (math.sqrt(((b)**2)-4.0*a*(c))))/(2.0*a)

print(k1)
print(k2)

A1=np.array([k1,-3*k1])
B1=np.array([5,k1])
C1=np.array([-k1,2])

A2=np.array([k2,-3*k2])
B2=np.array([5,k2])
C2=np.array([-k2,2])

dvec=np.array([-1,1])
omat=np.array([[0,1],[-1,0]])

A1B1=np.vstack((A1,B1)).T
B1C1=np.vstack((B1,C1)).T
C1A1=np.vstack((C1,A1)).T

A2B2=np.vstack((A2,B2)).T
B2C2=np.vstack((B2,C2)).T
C2A2=np.vstack((C2,A2)).T

def dir_vec(A1B1): return np.matmul(A1B1,dvec)
def norm_vec(A1B1): return np.matmul(omat,np.matmul(A1B1,dvec))

def dir_vec(A2B2): return np.matmul(A2B2,dvec)
def norm_vec(A2B2): return np.matmul(omat,np.matmul(A2B2,dvec))

def alt_point1(A1,B1C1):
	n1=dir_vec(B1C1)
	n2=norm_vec(B1C1)
	N=np.vstack((n1,n2))
	p=np.zeros(2)
	p[0]=np.matmul(n1,A1.T)
	p[1]=np.matmul(n2,B1C1[:,0])
	return np.matmul(np.linalg.inv(N),p)

def alt_point2(A2,B2C2):
	n1=dir_vec(B2C2)
	n2=norm_vec(B2C2)
	N=np.vstack((n1,n2))
	p=np.zeros(2)
	p[0]=np.matmul(n1,A2.T)
	p[1]=np.matmul(n2,B2C2[:,0])
	return np.matmul(np.linalg.inv(N),p)

def H1(A1,B1):
	n1=dir_vec(B1C1)
	n2=dir_vec(C1A1)
	N=np.vstack((n1,n2))
	p=np.zeros(2)
	p[0]=np.matmul(n1,A1.T)
	p[1]=np.matmul(n2,B1.T)
	return np.matmul(np.linalg.inv(N),p)

def H2(A2,B2):
	n1=dir_vec(B2C2)
	n2=dir_vec(C2A2)
	N=np.vstack((n1,n2))
	p=np.zeros(2)
	p[0]=np.matmul(n1,A2.T)
	p[1]=np.matmul(n2,B2.T)
	return np.matmul(np.linalg.inv(N),p)

P1=alt_point1(A1,B1C1)
Q1=alt_point1(B1,C1A1)
R1=alt_point1(C1,A1B1)
H1=H1(A1,B1)

P2=alt_point2(A2,B2C2)
Q2=alt_point2(B2,C2A2)
R2=alt_point2(C2,A2B2)
H2=H2(A2,B2)

print(H1)
print(H2)

len=10

lam_1 = np.linspace(0,1,len)
x_A1B1 = np.zeros((2,len))
x_B1C1 = np.zeros((2,len))
x_C1A1 = np.zeros((2,len))
x_A1P1 = np.zeros((2,len))
x_B1Q1 = np.zeros((2,len))
x_C1R1 = np.zeros((2,len))

x_A2B2 = np.zeros((2,len))
x_B2C2 = np.zeros((2,len))
x_C2A2 = np.zeros((2,len))
x_A2P2 = np.zeros((2,len))
x_B2Q2 = np.zeros((2,len))
x_C2R2 = np.zeros((2,len))

for i in range(len):
	temp1 = A1 + lam_1[i]*(B1-A1)
	x_A1B1[:,i]= temp1.T
	temp2 = B1 + lam_1[i]*(C1-B1)
	x_B1C1[:,i]= temp2.T
	temp3 = C1 + lam_1[i]*(A1-C1)
	x_C1A1[:,i]= temp3.T
	temp4 = A1 + lam_1[i]*(P1-A1)
	x_A1P1[:,i]= temp4.T
	temp5 = B1 + lam_1[i]*(Q1-B1)
	x_B1Q1[:,i]= temp5.T
	temp6 = C1 + lam_1[i]*(R1-C1)
	x_C1R1[:,i]= temp6.T



for i in range(len):
	tem1 = A2 + lam_1[i]*(B2-A2)
	x_A2B2[:,i]= tem1.T
	tem2 = B2 + lam_1[i]*(C2-B2)
	x_B2C2[:,i]= tem2.T
	tem3 = C2 + lam_1[i]*(A2-C2)
	x_C2A2[:,i]= tem3.T
	tem4 = A2 + lam_1[i]*5*(P2-A2)
	x_A2P2[:,i]= tem4.T
	tem5 = B2 + lam_1[i]*10*(Q2-B2)
	x_B2Q2[:,i]= tem5.T
	tem6 = R2 + lam_1[i]*14*(C2-R2)
	x_C2R2[:,i]= tem6.T

plt.subplot(2,1,1)
plt.plot(x_A1B1[0,:],x_A1B1[1,:],label='$A1B1$')
plt.plot(x_B1C1[0,:],x_B1C1[1,:],label='$B1C1$')
plt.plot(x_C1A1[0,:],x_C1A1[1,:],label='$C1A1$')
plt.plot(x_A1P1[0,:],x_A1P1[1,:],label='$A1P1$')
plt.plot(x_B1Q1[0,:],x_B1Q1[1,:],label='$B1Q1$')
plt.plot(x_C1R1[0,:],x_C1R1[1,:],label='$C1R1$')
plt.plot(A1[0], A1[1], 'o')
plt.text(A1[0] * (1 + 0.1), A1[1] * (1 + 0.05) , 'A1')
plt.plot(B1[0], B1[1], 'o')
plt.text(B1[0] * (1 + 0.02), B1[1] * (1) , 'B1')
plt.plot(C1[0], C1[1], 'o')
plt.text(C1[0] * (1 + 0.03), C1[1] * (1 - 0.1) ,'C1')
plt.plot(P1[0], P1[1], 'o')
plt.text(P1[0] * (1 + 0.1), P1[1] * (1 - 0.1) , 'P1')
plt.plot(Q1[0], Q1[1], 'o')
plt.text(Q1[0] * (1 + 0.1), Q1[1] * (1 + 0.4) , 'Q1')
plt.plot(R1[0], R1[1], 'o')
plt.text(R1[0] * (1 + 0.02), R1[1] * (1 + 0.2) , 'R1')
plt.plot(H1[0], H1[1], 'o')
plt.text(H1[0] * (1 + 0.1), H1[1] * (1 - 0.9) , 'H1')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid()

plt.subplot(2,1,2)
plt.plot(x_A2B2[0,:],x_A2B2[1,:],label='$A2B2$')
plt.plot(x_B2C2[0,:],x_B2C2[1,:],label='$B2C2$')
plt.plot(x_C2A2[0,:],x_C2A2[1,:],label='$C2A2$')
plt.plot(x_A2P2[0,:],x_A2P2[1,:],label='$A2P2$')
plt.plot(x_B2Q2[0,:],x_B2Q2[1,:],label='$B2Q2$')
plt.plot(x_C2R2[0,:],x_C2R2[1,:],label='$C2R2$')
plt.plot(A2[0], A2[1], 'o')
plt.text(A2[0] * (1 + 0.1), A2[1] * (1 - 0.1) , 'A2')
plt.plot(B2[0], B2[1], 'o')
plt.text(B2[0] * (1 - 0.2), B2[1] * (1) , 'B2')
plt.plot(C2[0], C2[1], 'o')
plt.text(C2[0] * (1 + 0.03), C2[1] * (1 - 0.9) ,'C2')
plt.plot(P2[0], P2[1], 'o')
plt.text(P2[0] * (1 + 0.1), P2[1] * (1 - 0.1) , 'P2')
plt.plot(Q2[0], Q2[1], 'o')
plt.text(Q2[0] * (1 + 0.1), Q2[1] * (1 - 0.9) , 'Q2')
plt.plot(R2[0], R2[1], 'o')
plt.text(R2[0] * (1 + 0.2), R2[1] * (1 - 0.8) , 'R2')
plt.plot(H2[0], H2[1], 'o')
plt.text(H2[0] * (1 + 0.02), H2[1] * (1 - 0.1) , 'H2')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid()
plt.show()
