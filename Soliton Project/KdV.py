#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 
from IPython.display import HTML


# In[2]:


class KdV(object):
    
    def __init__(self, dx = None, dt = None, eps = None, mu = None, x1 = None, xN = None, 
                 t1 = None, tN = None, anim_file_name = None, data_file_name = None):
        self.dx = dx
        self.dt = dt
        self.eps = eps    
        self.mu = mu
        self.x1 = x1
        self.xN = xN
        self.t1 = t1    
        self.tN = tN
        self.x = np.arange(x1, xN + dx, dx)
        self.t = np.arange(t1, tN + dt, dt)
        self.anim_file_name = anim_file_name
        self.data_file_name = data_file_name
        return
    def stability(self):
        return (self.dt/self.dx)*(self.eps + 4*self.mu/(self.dx**2))
    
    def initial_value(self):
        #u_t0 = (1/2)*(1 - np.tanh((self.x-(self.xN+self.x1))/5))
        
        #Gaussian wave packet
        #u_t0 = np.exp(-1/2*((self.x-(self.xN+self.x1))/2)**2)
        
        #Two soliton waves
        u_t0_1 = np.exp(-1/2*((self.x-(10))/2)**2)
        u_t0_2 = np.exp(-1/2*((self.x-(-10))/2)**2)*2
        u_t0 = u_t0_1 + u_t0_2
        return u_t0
    
    def initialize_matrix(self):
	#tanh
        #u_t0 = self.initial_value()
        #u = np.zeros((self.x.size + 2, self.t.size))
        #u[1:self.x.size+1,0] = u_t0
        #u[0,:] = u[1,0] 
        #u[1,:] = u[1,0]
        #u[self.x.size,:] = u[self.x.size,0] 
        #u[self.x.size+1,:] = u[self.x.size,0] 
	
	# sin^2
        u_t0 = np.zeros(self.x.size)
        for i in np.arange(self.x.size//5,2*self.x.size//5+1):
            u_t0[i] = np.sin(np.pi*((i-self.x.size//5)/(2*self.x.size//5-self.x.size//5)))**2
        u = np.zeros((self.x.size + 2, self.t.size))
        u[1:self.x.size+1,0] = u_t0
        u[0,:] = u[1,0] 
        u[1,:] = u[1,0]
        u[self.x.size,:] = u[self.x.size,0] 
        u[self.x.size+1,:] = u[self.x.size,0]
        return u

    
    def finite_diff(self):
        u = self.initialize_matrix()
        for i in np.arange(1, self.x.size + 1):
            if i == 1 or i == self.x.size:
                continue
            else:
                u[i,1] = u[i,0] -(self.eps/6)*(self.dt/self.dx)*(u[i+1,0]+u[i,0]+u[i-1,0])*(u[i+1,0]-u[i-1,0]) -(self.mu/2)*(self.dt/(self.dx**3))*(u[i+2,0]+2*u[i-1,0]-2*u[i+1,0]-u[i-2,0])

        
        for j in np.arange(2, self.t.size):
            for i in np.arange(1, self.x.size + 1):
                if i == 1 or i == self.x.size:
                    continue
                else:
                    u[i,j] = u[i,j-2] -(self.eps/3)*(self.dt/self.dx)*(u[i+1,j-1]+u[i,j-1]+u[i-1,j-1])*(u[i+1,j-1]-u[i-1,j-1]) -(self.mu)*(self.dt/(self.dx**3))*(u[i+2,j-1]+2*u[i-1,j-1]-2*u[i+1,j-1]-u[i-2,j-1])
        return u[1:-1,:]
    
    def animate(self):
        u = self.finite_diff()
        figure, ax = plt.subplots()
        
        # Setting limits for x and y axis
        ax.set_xlim(1.1*self.x1, 1.1*self.xN)
        ax.set_ylim(-0.1+np.min(u), 1.1*np.max(u))
        ax.set_xlabel('x')
        ax.set_ylabel('Height')
 
        # Since plotting a single graph
        line,  = ax.plot(0, 0) 
 
        def animation_function(i):
   
 
            line.set_xdata(self.x)
            line.set_ydata(u[:,i])
            return line,
 
        anim = animation.FuncAnimation(figure,
                          func = animation_function,
                          frames = np.arange(0, self.t.size, 1), 
                          interval = 10)

        HTML(anim.to_html5_video())
        writervideo = animation.FFMpegWriter(fps=60) 
        anim.save(self.anim_file_name, writer=writervideo) 
        plt.close()
        
        np.savetxt(self.data_file_name, u)


# In[3]:


class KdV_advection(object):
    
    def __init__(self, dx = None, dt = None, eps = None, mu = 0, local_value = None, x1 = None, xN = None, 
                 t1 = None, tN = None, anim_file_name = None, data_file_name = None):
        self.dx = dx
        self.dt = dt
        self.eps = eps    
        self.mu = mu
        self.x1 = x1
        self.xN = xN
        self.t1 = t1    
        self.tN = tN
        self.x = np.arange(x1, xN + dx, dx)
        self.t = np.arange(t1, tN + dt, dt)
        self.anim_file_name = anim_file_name
        self.data_file_name = data_file_name
        self.local_value = local_value
        return
    
    def stability(self):
        return (self.dt/self.dx)*(self.eps + 4*self.mu/(self.dx**2))
    
    def initialize_matrix(self):
        u_t0 = np.zeros(self.x.size)
        for i in np.arange(self.x.size//5,2*self.x.size//5+1):
            u_t0[i] = np.sin(np.pi*((i-self.x.size//5)/(2*self.x.size//5-self.x.size//5)))**2
        u = np.zeros((self.x.size + 2, self.t.size))
        u[1:self.x.size+1,0] = u_t0
        u[0,:] = u[1,0] 
        u[1,:] = u[1,0]
        u[self.x.size,:] = u[self.x.size,0] 
        u[self.x.size+1,:] = u[self.x.size,0] 
        return u
    
    def finite_diff(self):
        u = self.initialize_matrix()
        for i in np.arange(1, self.x.size + 1):
            if i == 1 or i == self.x.size:
                continue
            else:
                u[i,1] = u[i,0] -(self.eps/2)*(self.dt/self.dx)*(self.local_value)*(u[i+1,0]-u[i-1,0]) -(self.mu/2)*(self.dt/(self.dx**3))*(u[i+2,0]+2*u[i-1,0]-2*u[i+1,0]-u[i-2,0])

        
        for j in np.arange(2, self.t.size):
            for i in np.arange(1, self.x.size + 1):
                if i == 1 or i == self.x.size:
                    continue
                else:
                    u[i,j] = u[i,j-2] -(self.eps)*(self.dt/self.dx)*(self.local_value)*(u[i+1,j-1]-u[i-1,j-1]) -(self.mu)*(self.dt/(self.dx**3))*(u[i+2,j-1]+2*u[i-1,j-1]-2*u[i+1,j-1]-u[i-2,j-1])
        return u[1:-1,:]
    
    def animate(self):
        u = self.finite_diff()
        figure, ax = plt.subplots()
        
        # Setting limits for x and y axis
        ax.set_xlim(1.1*self.x1, 1.1*self.xN)
        ax.set_ylim(-0.1+np.min(u), 1.1*np.max(u))
        ax.set_xlabel('x')
        ax.set_ylabel('Height')
 
        # Since plotting a single graph
        line,  = ax.plot(0, 0) 
 
        def animation_function(i):
   
 
            line.set_xdata(self.x)
            line.set_ydata(u[:,i])
            return line,
 
        anim = animation.FuncAnimation(figure,
                          func = animation_function,
                          frames = np.arange(0, self.t.size, 1), 
                          interval = 10)

        HTML(anim.to_html5_video())
        writervideo = animation.FFMpegWriter(fps=60) 
        anim.save(self.anim_file_name, writer=writervideo) 
        plt.close()
        
        np.savetxt(self.data_file_name, u)


# In[ ]:




