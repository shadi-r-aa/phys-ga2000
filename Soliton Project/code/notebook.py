#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from KdV import KdV, KdV_advection


# In[2]:


kdv = KdV(dx = 0.4, dt = 0.05, eps = 0.2, mu = 0.1, x1 = -50, xN = 50, t1 = 0,
          tN = 200, anim_file_name = 'animation_tanh.mp4', data_file_name = 'data_tanh.txt')
kdv.stability()


# In[3]:


u = kdv.finite_diff()


# In[4]:


fig, axs = plt.subplots(5,1, sharex=True)
ax0, ax1, ax2, ax3, ax4 = axs[0], axs[1], axs[2], axs[3], axs[4]
fig.text(0.5, 0.04, 'x', ha='center')
fig.text(0.04, 0.5, 'U(x,t)', va='center', rotation='vertical')

    
ax0.plot(kdv.x,u[:,0],label = 't = 1')
ax0.legend()
ax1.plot(kdv.x,u[:,999],label = 't = 1000')
ax1.legend()
ax2.plot(kdv.x,u[:,1999],label = 't = 2000')
ax2.legend()
ax3.plot(kdv.x,u[:,2999],label = 't = 3000')
ax3.legend()
ax4.plot(kdv.x,u[:,3999],label = 't = 4000')
ax4.legend()
plt.savefig('tanh.png')


# In[5]:


kdv_advection = KdV_advection(dx = 0.4, dt = 0.05, eps = 0.2, local_value = 1, x1 = -50, xN = 50, t1 = 0,
          tN = 400, anim_file_name = 'animation_advection.mp4', data_file_name = 'data_advection.txt')
kdv_advection.stability()


# In[6]:


u = kdv_advection.finite_diff()


# In[7]:


fig, axs = plt.subplots(3,1, sharex=True)
ax0, ax1, ax2 = axs[0], axs[1], axs[2]
fig.text(0.5, 0.04, 'x', ha='center')
fig.text(0.04, 0.5, 'U(x,t)', va='center', rotation='vertical')

    
ax0.plot(kdv_advection.x,u[:,0],label = 't = 1')
ax0.legend()
ax1.plot(kdv_advection.x,u[:,1999],label = 't = 2000')
ax1.legend()
ax2.plot(kdv_advection.x,u[:,3999],label = 't = 4000')
ax2.legend()
plt.savefig('adv.png')


# In[ ]:




