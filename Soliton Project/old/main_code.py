#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt


# In[6]:


x_N = 100
x_1 = -100
t_N = 10
t_1 = 0

eps = 6
mu = 1

dx = 1
dt = 0.01

x = np.arange(x_1,x_N + 0.1, dx)
t = np.arange(t_1,t_N + 0.1, dt)

u_t0 = np.cos((np.pi*x)/(x_1*2))
plt.xlabel('Position')
plt.ylabel('Disturbance u')
plt.plot(x,u_t0)

u = np.zeros((x.size + 2, t.size))


u[1:x.size+1,0] = u_t0
u[0,:] = u[1,0] 
u[1,:] = u[1,0]
u[x.size,:] = u[x.size,0] 
u[x.size+1,:] = u[x.size,0] 


for i in np.arange(1, x.size + 1):
    if i == 1 or i == x.size:
        continue
    else:
        u[i,1] = u[i,0] -(eps/6)*(dt/dx)*(u[i+1,0]+u[i,0]+u[i-1,0])*(u[i+1,0]-u[i-1,0]) -(mu/2)*(dt/(dx**3))*(u[i+2,0]+2*u[i-1,0]-2*u[i+1,0]-u[i-2,0])

        
for j in np.arange(2, t.size):
    for i in np.arange(1, x.size + 1):
        if i == 1 or i == x.size:
            continue
        else:
            u[i,j] = u[i,j-2] -(eps/3)*(dt/dx)*(u[i+1,j-1]+u[i,j-1]+u[i-1,j-1])*(u[i+1,j-1]-u[i-1,j-1]) -(mu)*(dt/(dx**3))*(u[i+2,j-1]+2*u[i-1,j-1]-2*u[i+1,j-1]-u[i-2,j-1])


X, T = np.meshgrid(x, t)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, u[1:x.size+1, :].T, cmap='OrRd')  # Transpose u for correct dimensions
ax.set_xlabel('Position')
ax.set_ylabel('Time')
ax.set_zlabel('Disturbance u')
plt.show()


# In[ ]:




