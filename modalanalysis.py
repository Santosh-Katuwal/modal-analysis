import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
k=[1800,1200,600,600,600,600,500,300,300,300,300,300,300,300,300] #enter floor stiffness from bottom to top [k1,k2,k3...]
m=[2.5,2,1.5,1,1,1,1,1,1,1,1,1,1,1,1]         #enter floor mass from bottom to top [m1,m2,m3...]
trials=20        #Define number of trials to be carriedout
n=len(k)        #nos of DOF
K0=np.zeros((n,n))
M0=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i==j and i<n-1:
            try:
                K0[i,j]=k[i]+k[i+1]
                M0[i,j]=m[i]
            except:
                pass
        elif i==j and i==n-1:
            K0[i,j]=k[i]
            M0[i,j]=m[i]
        elif abs(i-j)==1 and i>j:
            K0[i,j]=-k[i]
        elif abs(i-j)==1 and i<n:
            K0[i,j]=-k[j]
        else:
            pass
#Flexibility matrix
f=np.linalg.inv(K0)
D=f*m
phi_1=np.ones((n,1))    #Trial mode shape

DM=[]
PHI=[]
OMGEGA_i=[]
for i in range(1,n+1):  #i =mode number
    if i==1:
       for j in range(0,trials):
           Phi_next=np.matmul(D,phi_1)
           phi_next=Phi_next/np.abs(Phi_next).max()
           phi_1=phi_next
       Dphi=np.matmul(D,phi_1)
       omega_1=np.sqrt(np.abs(phi_1).max()/np.abs(Dphi).max())
       X=np.matmul(phi_1,np.matmul(phi_1.transpose(),M0))  #Assume X=phi_1 * phi_1^T  * mass
       Y=np.matmul(phi_1.transpose(),np.matmul(M0,phi_1))  #Assume Y=phi_1^T *  mass *phi_1 
       I=np.identity(n)
       S1=I-X/Y
       Si=S1
    if i>1:
        dm=globals()[f"D{i}"]=np.matmul(D,Si)
        DM.append(dm)
        del dm #Deleting variable after appending to reduce memory use
        phi_new=np.ones((n,1))    #Trial mode shape
        ph=globals()[f"phi_{i}"]=phi_new
        PHI.append(ph)
        del ph
        for j in range(0,trials):
            Phi_next=np.matmul(DM[i-2],PHI[i-2])   
            phi_next=Phi_next/np.abs(Phi_next).max() 
            phi_new=phi_next
            PHI[i-2]=phi_new
        Dphi_i=np.matmul(DM[i-2],PHI[i-2])
        omega_i=np.sqrt(np.abs(PHI[i-2]).max()/np.abs(Dphi_i).max())
        OMGEGA_i.append(omega_i)
        X=np.matmul(PHI[i-2],np.matmul(PHI[i-2].transpose(),M0))  #Assume X=phi_1 * phi_1^T  * mass
        Y=np.matmul(PHI[i-2].transpose(),np.matmul(M0,PHI[i-2]))  #Assume Y=phi_1^T *  mass *phi_1 
        Si=Si-X/Y

floors=[0]
flrs=np.linspace(0, 3, 10)
h=3
for i in range(len(k)):
    floors.append(h+h*i)


#printing
disp_phi1=(phi_1.round(3))
new_phi1=np.insert(disp_phi1,0,0)
cir_fr1=omega_1



disp_phi2=(PHI[0].round(3))
new_phi2=np.insert(disp_phi2,0,0)
cir_fr2=OMGEGA_i[0]


disp_phi3=(PHI[1].round(3))
new_phi3=np.insert(disp_phi3,0,0)
cir_fr3=OMGEGA_i[1]


disp_phi4=(PHI[2].round(3))
new_phi4=np.insert(disp_phi4,0,0)
cir_fr4=OMGEGA_i[2]

disp_phi5=(PHI[3].round(3))
new_phi5=np.insert(disp_phi5,0,0)
cir_fr5=OMGEGA_i[3]



fig,axs=plt.subplots(1,3)

loops=np.arange(1,50,0.01)
for t in loops:
    axs[0].cla()
    A=new_phi1*np.sin(cir_fr1*t)
    axs[0].scatter(A-0.5,floors,s=50,c='r')
    axs[0].scatter(A+0.5,floors,s=50,c='r')
    for i in range(len(k)):
        axs[0].plot([A[i]-0.5,A[i+1]-0.5],[floors[i],floors[i+1]],c='k')
        axs[0].plot([A[i+1]-0.5,A[i+1]+0.5],[floors[i+1],floors[i+1]],c='k')
        axs[0].plot([A[i]+0.5,A[i+1]+0.5],[floors[i],floors[i+1]],c='k')
    axs[0].set_xlim(-2,2)
    axs[0].set_ylim(0,max(floors)+1)
    plt.tight_layout()
    #plt.pause(0.01)

    axs[1].cla()
    A2=new_phi2*np.sin(cir_fr2*t)
    axs[1].scatter(A2-0.5,floors,s=50,c='r')
    axs[1].scatter(A2+0.5,floors,s=50,c='r')
    for i in range(len(k)):
        axs[1].plot([A2[i]-0.5,A2[i+1]-0.5],[floors[i],floors[i+1]],c='k')
        axs[1].plot([A2[i+1]-0.5,A2[i+1]+0.5],[floors[i+1],floors[i+1]],c='k')
        axs[1].plot([A2[i]+0.5,A2[i+1]+0.5],[floors[i],floors[i+1]],c='k')
    axs[1].set_xlim(-2,2)
    axs[1].set_ylim(0,max(floors)+1)
    plt.tight_layout()

    axs[1].scatter(A2-0.5,floors,s=50,c='r')
    plt.scatter(A2+0.5,floors,s=50,c='r')
    for i in range(len(k)):
        axs[1].plot([A2[i]-0.5,A2[i+1]-0.5],[floors[i],floors[i+1]],c='k')
        axs[1].plot([A2[i+1]-0.5,A2[i+1]+0.5],[floors[i+1],floors[i+1]],c='k')
        axs[1].plot([A2[i]+0.5,A2[i+1]+0.5],[floors[i],floors[i+1]],c='k')
    axs[1].set_xlim(-2,2)
    axs[1].set_ylim(0,max(floors)+1)
    plt.tight_layout()
    
    axs[2].cla()
    A3=new_phi3*np.sin(cir_fr3*t)
    axs[2].scatter(A3-0.5,floors,s=50,c='r')
    plt.scatter(A3+0.5,floors,s=50,c='r')
    for i in range(len(k)):
        axs[2].plot([A3[i]-0.5,A3[i+1]-0.5],[floors[i],floors[i+1]],c='k')
        axs[2].plot([A3[i+1]-0.5,A3[i+1]+0.5],[floors[i+1],floors[i+1]],c='k')
        axs[2].plot([A3[i]+0.5,A3[i+1]+0.5],[floors[i],floors[i+1]],c='k')
    axs[2].set_xlim(-2,2)
    axs[2].set_ylim(0,max(floors)+1)
    plt.tight_layout()

    
    plt.pause(0.01)
plt.xlabel('Relative Amplitude')
plt.ylabel('Floor height')

