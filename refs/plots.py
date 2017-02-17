
# coding: utf-8

# In[1]:

import scipy as sp
import scipy.constants
import scipy.integrate
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import numpy as np
import time


# In[94]:

data = np.loadtxt("../data/test.txt")


# In[95]:

plt.plot(data[:, 0], data[:, 1], "-")
plt.axhline(0.0, linestyle="--", color="k")
#plt.axvline(1.64692)
#plt.xlim([0.99, 1.01])
#plt.ylim([0.0,1.0])
plt.show()


# In[ ]:




# In[75]:

0.0 5.86362e-047
-1.0 5.70575e-047
-2.0 5.59936e-047
-3.0 5.62015e-047
-4.0 5.62985e-047
-5.0 5.64122e-047
-6.0 5.62728e-047
-7.0 5.62192e-047
...


# In[102]:

data = np.loadtxt("../data/outf.txt")


# In[104]:

plt.plot(data[:, 0], data[:, 1], "-")
#plt.axvline(1.64692)
#plt.xlim([0.99, 1.01])
#plt.ylim([0.0,1.0])
#plt.yscale('log')
plt.show()


# In[ ]:



