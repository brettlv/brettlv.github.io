import numpy as np
import random
import math
from astropy.constants import kpc,c,pc
import matplotlib.pyplot as plt

# def mean():
#     data = np.loadtxt("data.txt")
#     mjd = data[:,0]
#     "{},{},{}\n".format(i,declination,ascension))
#     newfile.closenewfile = open("mean.txt",'w')
#     for i in set (mjd):
#         s = np.where(mjd==i)
#         declination = np.mean(data[s[0],1])
#         ascension = np.mean(data[s[0],2])
#         newfile.write("{},{},{}\n".format(i,declination,ascension))
#     newfile.close

def read_file():
    mjd_mean=[]
    de_obs=[]
    as_obs=[]
    mjd_mean,de_obs,as_obs = np.loadtxt("data.txt",unpack=True)
    return mjd_mean,de_obs,as_obs

def fun(varable):
    mjd_mean,as_obs,de_obs = read_file()
    aim = 0
    L = 1.907 * pc.value *10**9
    c_day = c.value*24*60*60
    i = varable[0]
    R_A = varable[1]
    lamda = varable[2]
    epsilon = varable[3]
    omega = varable[4]
    eta = varable[5]
    for j in range(len(mjd_mean)):
        t0 = mjd_mean[j]-mjd_mean[0]
        eta_a = eta + omega*t0
        R_y = R_A * (np.sin(lamda) * np.sin(eta_a))
        R_z = R_A * (np.cos(lamda) * np.sin(i) - np.sin(lamda) * np.cos(i) * np.cos(eta_a))
        Ra = (np.degrees((R_y * np.sin(epsilon) + R_z * np.cos(epsilon))/L))*60*60*10**3
        Rb = (np.degrees((R_y * np.cos(epsilon) - R_z * np.sin(epsilon))/L))*60*60*10**3
        aim1=(Ra-as_obs[j])**2 + (Rb - de_obs[j])**2
        aim += aim1
    return aim/4


def randomnum(n):
    x=np.random.choice(n,1)
    return x

def deltaf(y,k):
    delta = 0
    c = 1000
    if y==0:
        if k>c:
            delta = 0.0001
        else:
            delta = 0.001
    elif y==1:
        if k >c:
            delta = 1 * pc.value
        else:
            delta = 0.1 * pc.value
    elif y==2:
        if k>c:
            delta = 0.0001
        else:
            delta = 0.001
    elif y==3:
        if k>c:
            delta = 0.01
        else:
            delta = 0.1
    elif y==4:
        if k>c:
            delta = 0.001
        else:
            delta = 0.005
    elif y==5:
        if k>c:
            delta = 0.01
        else:
            delta = 0.1
    return delta

def rangex(ptrial,k,x,delta,y):
    while True:
        ptrial[k,y] = x[k,y]+delta*(2.0*random.random()-1)
        if y == 0 and (0.008<ptrial[k,0]<0.015):
            break
        elif y == 1 and (320*pc.value<ptrial[k,1]<330*pc.value):
            break
        elif y == 2 and (0.006<ptrial[k,2]<0.015):
            break
        elif y == 3 and (4.5<ptrial[k,3]<5.5):
            break
        elif y == 4 and (0.008<ptrial[k,4]<0.015):
            break
        elif y == 5 and (0.9<ptrial[k,5]<1.5):
            break
    return

def initial(x):
    r=random.random()
    x[0,0] = 0.008+0.001*r
    x[0,1] = 320.0*pc.value+1*pc.value*r
    x[0,2] = 0.008+0.001*r
    x[0,3] = 4.95 + 0.1*r
    x[0,4] = 0.01 + 0.005*r
    x[0,5] = 1.0 + 0.01*r
    return

def SA(ptrial,n,k,fc, deltaE_avg,t, number_accept):

    for j in range(n):
        fcn[j] = ptrial[k,j]
    f_new = fun(fcn)
    f_old = fun(fc)
    for j in range(n):
        best_v[j] = fcn[j]
    f_best = min(f_old,f_new)
    deltaE = abs(f_old-f_new)
    if k == 0:
        deltaE_avg = deltaE
    ##比较得出最优解
    if f_new > f_best:
        p = math.exp(-deltaE/(deltaE_avg * t))
        print(deltaE, deltaE_avg ,t, p)
        if (random.random() < p):
            accept = True
        else:
            accept = False
    else:
        accept = True
    if accept :
        f_best = f_new
        number_accept = number_accept + 1.0
        deltaE_avg = (deltaE_avg * (number_accept - 1.0) + deltaE) / number_accept
        for j in range(n):
            best_v[j] = fcn[j]
    else:
        for j in range(n):
            best_v[j] = fc[j]
    t = frac * t
    for j in range(n):
        ptrial[k+1,j] = best_v[j]
        x[k+1,j] = best_v[j]
    return (f_best,best_v,deltaE_avg, number_accept)

def calculate_Ra_Rb(best_v):
    L = 1.907 *10 **9 * pc.value
    i = best_v[0]
    R_A = best_v[1]
    lamda = best_v[2]
    epsilon = best_v[3]
    omega = best_v[4]
    eta = best_v[5]
    eta_a = eta
    R_y = R_A * (np.sin(lamda) * np.sin(eta_a))
    R_z = R_A * (np.cos(lamda) * np.sin(i) - np.sin(lamda) * np.cos(i) * np.cos(eta_a))
    Ra = (np.degrees((R_y * np.sin(epsilon) + R_z * np.cos(epsilon))/L))*60*60*10**3
    Rb = (np.degrees((R_y * np.cos(epsilon) - R_z * np.sin(epsilon))/L))*60*60*10**3
    return Ra, Rb

def draw(best_v):
    mjd=[]
    as_obs1=[]
    de_obs1=[]
    mjd,as_obs1,de_obs1=read_file()
    t0 = np.linspace(0,2*np.pi/best_v[4],1000)
    L = 1.907 * 10 **9 * pc.value
    i = best_v[0]
    R_A = best_v[1]
    lamda = best_v[2]
    epsilon = best_v[3]
    omega = best_v[4]
    eta = best_v[5]
    eta_a = eta+t0*omega
    R_y = R_A * (np.sin(lamda) * np.sin(eta_a))
    R_z = R_A * (np.cos(lamda) * np.sin(i) - np.sin(lamda) * np.cos(i) * np.cos(eta_a))
    Ra = (np.degrees((R_y * np.sin(epsilon) + R_z * np.cos(epsilon))/L))*60*60*10**3
    Rb = (np.degrees((R_y * np.cos(epsilon) - R_z * np.sin(epsilon))/L))*60*60*10**3
    plt.figure(figsize=(8,8))
    plt.plot(Ra,Rb,"b--")
    plt.scatter(as_obs1,de_obs1,c='r')
    plt.xlabel("Relative right Ascension")
    plt.ylabel("Relative declination")
    plt.title("Fitted trajectories of component of PKS 1510-089")
    plt.show()

##m为循环的次数,n为变量的个数, t0为初始温度.
delta = 10**(-3)
m = 1000000
n = 6
x = np.zeros((m+1,n),np.float)
ptrial = np.empty((m+1,n),np.float)
fc = np.empty((n),np.float)
fcn = np.empty((n),np.float)
best_v = np.empty((n),np.float)
deltaE_avg = 0
p1 = 0.7
pn = 0.01
t1 = -1/math.log(p1)
tn = -1/math.log(pn)
t = -1/np.log(0.7)
frac = (tn/t1)**(1.0/(m-1.0))
initial(x)
number_accept = 0
for k in range(m):
    for j in range(n):
        fc[j]=x[k,j]
        ptrial[k,j]=x[k,j]
    y=randomnum(n)
    delta=deltaf(y,k)
    rangex(ptrial,k,x,delta,y)
    f_best,best_v,deltaE_avg, number_accept = SA(ptrial,n,k,fc,deltaE_avg, t, number_accept)
    t = t*frac
    if f_best < 0.009:
        break
    print('Best solution: ',f_best)
    #print ('Best parameter: ',best_v)
print('Best',best_v)
print("calculate",calculate_Ra_Rb(best_v))
#### draw picture
draw(best_v)
