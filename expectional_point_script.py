# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 20:23:53 2022

@author: CavaSCXRD
"""

'''
Import Packages
'''
from qutip import *
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import cmath as cm

###########################################################################

'''
Experimental Set up
'''
#experimental parameters
w = 5.8*2*np.pi*10**3
gamma_e = 5.25
gamma_f = 0.25

#define the energy levels
g = basis(3,0)
e = basis(3,1)
f = basis(3,2)

#pauli matrices
sx = tensor(e,f.dag())+tensor(f,e.dag())
sy = -complex(0,1)*tensor(e,f.dag()) + complex(0,1)*tensor(f,e.dag())
sz = tensor(e,e.dag()) - tensor(f,f.dag())

'''
#effective hamiltonian
H_eff = J*(tensor(e,f.dag())+tensor(f,e.dag()))-(delta/2)*(tensor(f,f.dag())-tensor(e,e.dag()))
print(H_eff)
H_eff.dims = [[3],[3]]
'''

#initial state
psi0 = basis(3,2) #start in f state
#print(psi0)

#Dissipation Parameters
#rate from e-state
L_e = np.sqrt(gamma_e)*tensor(g,e.dag())
L_e.dims = [[3],[3]]

#rate from f-state
L_f = np.sqrt(gamma_f)*tensor(e,f.dag())
L_f.dims = [[3],[3]]

d_ops = [L_e, L_f]

###########################################################################

'''
Simulation
'''

def Evolution(psi0, J=[0], delta=[0], time=2, dissipation=d_ops):
    #Definitions
    #psi0 is the initial state
    #J is the coupling
    #delta is the detuning
    #time is the simulation time in microseconds
    #dissipation is the 'emmision' array
    #----------------------------------------------
    
    #how long we're running simulation
    tlist = np.linspace(0,time,1000) #2 microseconds
    
    #initialize all data array
    all_res=[]
    
    for i in delta:
        for j in J:
            
            #effective hamiltonian
            H_eff = j*(tensor(e,f.dag())+tensor(f,e.dag()))-(i/2)*(tensor(f,f.dag())-tensor(e,e.dag()))
            #H_eff = j*(tensor(e,f.dag())+tensor(f,e.dag()))+(i-(complex(0,1)*gamma_e)/2)*tensor(e,e.dag())
            #print(H_eff)
            H_eff.dims = [[3],[3]]
            
            #Population result
            result = mesolve(H_eff, psi0, tlist, dissipation)
            
            #raw data
            rho_ee = [state[1][0][1].real for state in result.states]
            rho_ff = [state[2][0][2].real for state in result.states]
            
            #normalized data
            P_ee_n = np.array(rho_ee)/(np.array(rho_ff)+np.array(rho_ee))
            P_ff_n = np.array(rho_ff)/(np.array(rho_ff)+np.array(rho_ee))
            
            #cartesian result
            
            
            #Pauli matrices act on bloch sphere, but we expand the matrices to 
            sx.dims=[[3],[3]]
            sy.dims=[[3],[3]]
            sz.dims=[[3],[3]]
            
            
            cart_result = mesolve(H_eff, psi0, tlist, dissipation, [sx, sy, sz])
            
            
            
            data = [psi0, i, j, P_ee_n, P_ff_n, cart_result.expect[0], cart_result.expect[1], cart_result.expect[2]]
            all_res.append(data)
    
    return all_res





###########################################################################

'''
Theory
'''

def Theory(time=2, J=1):
    
    #how long we're running simulation
    tlist = np.linspace(0,time,1000) #2 microseconds
    
    alpha = np.sqrt(J*J - gamma_e**2/16)
    theta = np.arcsin(gamma_e/4/J)
    
    #evolution of Pe in time
    def Pe(time):
        return np.exp(-gamma_e/2*time)*(J/alpha)**2*np.sin(alpha*time)**2

    #evolution of Pf in time
    def Pf(time):
        return np.exp(-gamma_e/2*time)*(J/alpha)**2*np.cos(alpha*time-theta)**2
    
    pop_e=[]
    pop_f=[]
    for time in tlist:
        eee=Pe(time)
        fff=Pf(time)
        tot=eee+fff
        #append normlaized values
        pop_e.append(eee/tot)
        pop_f.append(fff/tot)
    
    
    return [pop_e, pop_f]





###########################################################################

'''
Plotting section
'''

def Plot_data(image, x=[], y=[], z=[], ax=[r'$t\;[\mu s]$', r'$P_{f}^{n}(t)$ $[rad\;\mu\;s^{-1}]$', 'z'], lbls=[]):
    if lbls==[] and len(x)>1:
        lbls=len(x)*[0]
    if image=='standard':
        fig, axes = plt.subplots(1,1)
        for i in range(len(x)):
            axes.plot(x[i], y[i], label=lbls[i])
        axes.set_xlabel(ax[0], fontsize=20)
        axes.set_ylabel(ax[1], fontsize=16)
    elif image == 'cmap':
        fig, axes = plt.subplots(1,1)
        fig.set_size_inches(6, 5)
        cmap = axes.pcolormesh(x, y, z)
        #axes.axhline(xmin=0.0,xmax=20.0,y=γ_e/4,linestyle='--',color='black', label = r'$EP=\gamma_{e}/4$')
        plt.gca().invert_yaxis()
        fig.colorbar(cmap)
        axes.set_xlabel(r'$\theta$', fontsize=20)
        axes.set_ylabel(r"$J\;[rad\;\mu s^{-1}]$", fontsize=20);
    elif image == 'bloch':
        fig = Bloch()
        fig.zlabel = [ax[0], ax[1]]
        fig.clear()
        for i in range(len(x)):
            fig.add_points([x[i], y[i], z[i]])
        #fig.add_points([sy_s, sy_s, sz_s])
        #fig.add_points([sx_b, sy_b, sz_b])
        fig.show()
    if image!='bloch':
        fig.legend()
    plt.show(fig)
    return