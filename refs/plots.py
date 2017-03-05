
# coding: utf-8

# In[1]:

import scipy as sp
import scipy.constants
import scipy.integrate
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np
import time


# In[13]:

params = ["_f4.0_p4.5_m7.0", "_f4.0_p4.5_m0.0"]

tunnel_data_arr = []; pot_arr = []

for par in params:
    data_file = "../data/tunnel_data"+par+".txt"
    pot_file = "../data/potential"+par+".txt"

    tunnel_data_arr.append(np.loadtxt(data_file))
    pot_arr.append(np.loadtxt(pot_file))


# In[14]:

for i, pot in enumerate(pot_arr):
    plt.plot(pot[:, 0], pot[:, 1], "-")
plt.axhline(0.0, linestyle="--", color="k")
#plt.axvline(1.64692)
#plt.xlim([0.99, 1.01])
#plt.ylim([0.0,1.0])
plt.show()


# In[19]:

plt.plot(tunnel_data_arr[0][:, 0], tunnel_data_arr[0][:, 1], "-")
plt.plot(tunnel_data_arr[1][:, 0]+7.0, tunnel_data_arr[1][:, 1], "-")
#plt.axvline(1.64692)
#plt.xlim([6.0, 11.0])
#plt.ylim([1e-6,1.0])
#plt.yscale('log')
plt.show()


# In[ ]:



