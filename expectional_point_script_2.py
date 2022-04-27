# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 01:49:36 2022

@author: CavaSCXRD
"""

#copy and paste this into terminal to run the 3b
i=np.pi/2
phi=np.linspace(0,2*np.pi,15)
J=np.linspace(0,5,15)
phi_res=[]
for j in phi:
    res=Evolution(np.cos(i/2)*basis(3,1)-complex(0,1)*np.exp(complex(0,1)*j)*np.sin(i/2)*basis(3,2), J, time=0.5)
    phi_res.append(res)
z=[]
for row in phi_res:
    temp_z=[]
    for i in row:
        temp_z.append(i[4][0]-i[4][-1])
    z.append(temp_z)
Plot_data('cmap', phi, J, np.transpose(z))

#copy and paste this into terminal to run 3a
j=np.pi/2
theta=np.linspace(0,np.pi,15)
J=np.linspace(0,1.5,15)
theta_res=[]
for i in theta:
    res=Evolution(np.cos(i/2)*basis(3,1)-complex(0,1)*np.exp(complex(0,1)*j)*np.sin(i/2)*basis(3,2), J, time=0.5)
    theta_res.append(res)
y=[]
k=0
for row in theta_res:
    temp_y=[]
    for i in row:
        temp_y.append(abs(i[4][0]-i[4][-1])+(i[4][0]-i[4][-1])/2)
    y.append(temp_y)
    k=k+1
Plot_data('cmap', theta, J, np.transpose(y))