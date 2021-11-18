#!/usr/bin/env python
# coding: utf-8

# # Solusi UTS Komputasi Sains

# 1.

# In[1]:


import numpy as np


# In[2]:


A = np.array([[1, 0, -1],
              [-0.5, 1, -0.25],
              [1, -0.5, 1]])


# In[3]:


def diagonal(A):
  n = A.shape[0]
  D = np.zeros_like(A)
  for i in range(n):
    D[i,i] = A[i,i]
  return D

def lower(A):
  n = A.shape[0]
  L = np.zeros_like(A)
  for i in range(n):
    for j in range(n):
      if i > j:
        L[i,j] = -A[i,j]
  return L

def upper(A):
  n = A.shape[0]
  U = np.zeros_like(A)
  for i in range(n):
    for j in range(n):
      if i < j:
        U[i,j] = -A[i,j]
  return U


# In[4]:


def spectral_radius(T):
  eigval = np.linalg.eigvals(T)
  rho = np.max(np.abs(eigval))
  return rho


# In[5]:


D = diagonal(A)
L = lower(A)
U = upper(A)
print("{} \n\n {} \n\n {}".format(D, L, U))


# In[6]:


Tj = np.matmul(D, (L+U))
Tj


# In[7]:


spectral_radius(Tj)


# In[8]:


X = np.linalg.inv(D - L)
Tgs = np.matmul(X, U)
Tgs


# In[9]:


spectral_radius(Tgs)


# 3.

# In[10]:


from scipy.linalg import norm


# In[11]:


def jacobi1(A, b, x0, epsilon=1e-5, N=1000):
  x = np.zeros(len(b))
  T = A - np.diag(np.diagonal(A))

  for k in range(N):
    x = (b - np.dot(T,x0))/np.diagonal(A)

    if norm(np.dot(A,x) - b) < epsilon:
      break
      
    x0 = x

  return k,x


# In[12]:


def gauss_seidel(A, b, x0, epsilon=1e-5, N=1000):
  x = np.zeros_like(b, dtype=np.double)
  
  for k in range(N):

    for i in range(len(b)):
      U = np.dot(A[i,:i], x[:i])
      V = np.dot(A[i,(i+1):], x0[(i+1):])
      x[i] = 1/A[i,i] * (b[i] - U - V)
      #print(k,x)
      
    if norm(np.dot(A,x) - b) < epsilon:
      break
    
    x0 = x

  return k,x


# a.

# In[13]:


A = np.array([[2, -1, 1],
              [2, 2, 2],
              [-1, -1, 2]], dtype = np.double)
b = np.array([-1, 4, -5], dtype = np.double)


# In[14]:


HasilGS = gauss_seidel(A, b, x0=np.zeros(A.shape[0]))
HasilGS


# In[15]:


HasilJ = jacobi1(A, b, x0=np.zeros(A.shape[0]))
HasilJ


# Kenapa?

# In[16]:


D = diagonal(A)
L = lower(A)
U = upper(A)

print("{} \n\n {} \n\n {}".format(D, L, U))


# * Untuk Jacobi

# $T_J = D^{-1}(L + U)$

# In[17]:


D_invers = np.linalg.inv(D)
Tj = np.matmul(D, (L + U))
Tj


# In[18]:


spectral_radius(Tj)


# Karena $\rho (T_J) > 1$ maka menghasilkan barisan solusi yang __divergen__.

# * Untuk Gauss-Seidel

# $T_{GS} = (D - L)^{-1} U$

# In[19]:


DL_invers = np.linalg.inv((D-L))
Tgs = np.matmul(DL_invers, U)
Tgs


# In[20]:


spectral_radius(Tgs)


# Karena $\rho (T_{GS}) < 1$ maka menghasilkan barisan solusi yang __konvergen__.

# b.

# In[21]:


A = np.array([[1, 2, -2],
              [1, 1, 1],
              [2, 2, 1]], dtype=np.double)
b = np.array([7, 2, 5], dtype=np.double)


# In[22]:


HasilGS = gauss_seidel(A, b, x0=np.zeros(A.shape[0]))
HasilGS


# In[23]:


HasilJ = jacobi1(A, b, x0=np.zeros(A.shape[0]))
HasilJ


# Kenapa?

# In[24]:


D = diagonal(A)
L = lower(A)
U = upper(A)

print("{} \n\n {} \n\n {}".format(D, L, U))


# * Untuk Jacobi

# In[25]:


D_invers = np.linalg.inv(D)
Tj = np.matmul(D, (L + U))
Tj


# In[26]:


spectral_radius(Tj)


# Karena $\rho (T_J) < 1$ maka menghasilkan barisan solusi yang __konvergen__.

# * Untuk Gauss-Seidel

# In[27]:


DL_invers = np.linalg.inv((D-L))
Tgs = np.matmul(DL_invers, U)
Tgs


# In[28]:


spectral_radius(Tgs)


# Karena $\rho (T_{GS}) > 1$ maka menghasilkan barisan solusi yang __divergen__.

# 4.

# In[29]:


import matplotlib.pyplot as plt


# In[30]:


def cubic_spline(data):
  n = len(data) - 1

  data = np.array(data)

  # Nilai dari koef. a diketahui dari titik data atau input yaitu a_j = f(x_j)
  a = [data[i, 1] for i in range(n+1)]
  a = np.array(a)

  h = [(data[i+1,0] - data[i,0]) for i in range(n)]
  
  # Membentuk matriks A
  A = np.zeros((n+1,n+1))
  for i in range(1, n):
    for j in range(0, n):
      if j < i:
        for k in range(j):
          A[i,k] = 0
        A[i,j] = h[j-1]
      elif j > i:
        for k in range(j):
          A[k,j] = 0
        A[i,j] = h[j]
        A[i+1, j+1] = h[j]
      else:
        A[i,i] = 2*(h[j-1] + h[j])
  A[0, 0] = 1
  A[n, n] = 1
  A[:n-1, n] = 0

  # Membentuk vektor b
  b = np.zeros(n+1)
  for i in range(1,n):
    b[i] = (3/h[i]) * (a[i+1] - a[i]) - (3/h[i-1]) * (a[i] - a[i-1])

  # Cari nilai koef. c dengan Gauss-Seidel
  x0 = np.zeros(n+1)
  iter, c = gauss_seidel(A, b, x0)

  # Mencari koef. b dan d setelah koef. c didapatkan
  d = np.zeros(n+1)
  for i in range(n):
    b[i] = (a[i+1] - a[i])/h[i] - h[i]*(c[i+1] + 2*c[i])/3
    d[i] = (c[i+1] - c[i])/3*h[i]

  return a, b, c, d, A


# In[31]:


x = np.array([0.9, 1.3, 1.9, 2.1, 2.6, 3.0, 3.9, 4.4, 4.7, 5.0, 6.0, 7.0, 8.0, 9.2, 10.5, 11.3, 11.6, 12.0, 12.6, 13.0, 13.3])
y = np.array([1.3, 1.5, 1.85, 2.12, 2.6, 2.7, 2.4, 2.15, 2.05, 2.11, 2.25, 2.3, 2.26, 1.95, 1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25])

data = [list(a) for a in zip(x, y)]
data = np.array(data)
data


# In[32]:


# Cari koef. untuk S
import pprint
pp = pprint.PrettyPrinter(indent=2, compact=True)
a,b,c,d, A = cubic_spline(data)

P = [a,b,c,d]
pp.pprint(P)


# In[33]:


n = len(data)
S = [(lambda x, j = j: a[j] + b[j]*(x - data[j,0]) + c[j]*((x - data[j,0])**2) + d[j]*((x - data[j,0])**3)) for j in range(n-1)]


# In[34]:


interval = [np.linspace(data[i,0], data[i+1, 0], 100) for i in range(n-1)]
plt.scatter(data[:,0], data[:,1], c='red', label='titik data')
for i in range(n-1):
  plt.plot(interval[i], S[i](interval[i]), 'b')

plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.ylim([-4, 4]);
plt.legend();


# In[ ]:




