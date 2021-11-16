#!/usr/bin/env python
# coding: utf-8

# # Metode Iteratif untuk Menyelesaikan Sistem Persamaan Linear
# 
# Terdapat dua teknik untuk menyelesaikan Sistem Persamaan Linear (SPL) yaitu metode langsung dan metode iteratif. Metode langsung seperti eliminasi Gauss, memiliki beberapa kelemahan diantaranya:
# * Tidak memiliki perbaikan dari error yang dihasilkan. Hal ini dikarenakan, metode ini terlalu mirip dengan solusi eksaknya. Akibatnya, metode ini rentan atau sensitif terhadap error pembulatan. 
# * Operasi yang dibutuhkan oleh metode langsung dinotasikan dengan $O(n^3)$. Ini berarti, semakin besar $n$ atau ukuran matriksnya, maka akan semakin lama waktu eksekusinya.

# ## Metode Jacobi
# 
# Metode Jacobi dibangun dengan menggunakan dua asumsi:
# 1. SPL $A\vec{x}=\vec{b}$ memiliki solusi tunggal.
# 2. Koefisien dari matriks $A$ tidak nol pada diagonal utamanya, $a_{11}, a_{12}, \cdots, a_{nn}$. Jika koefisien dari diagonal utamanya bernilai nol, maka penukaran baris atau kolom dapat dilakukan.
# 
# Misalkan diberikan SPL
# 
# $$
# \begin{array}{cc}
#   a_{11}x_1 & + & a_{12}x_2 & + & \cdots & + & a_{1n}x_n & = & b_1 & (1)\\
#   a_{21}x_1 & + & a_{22}x_2 & + & \cdots & + & a_{2n}x_n & = & b_2 & (2)\\
#    &  &  &  & \vdots &  &  &  &  \\
#   a_{n1}x_1 & + & a_{n2}x_2 & + & \cdots & + & a_{nn}x_n & = & b_n & (n)\\
# \end{array}
# $$
# 
# Berikut ini adalah langkah-langkah penyelesaian menggunakan metode Jacobi:
# 1. Pertama, kita selesaikan persamaan $(1)$ untuk $x_1$, dilanjutkan ke persamaan $(2)$ untuk $x_2$, dan seterusnya sampai persamaan $(n)$ untuk $x_n$. Sehingga menjadi
# 
# $$
# \begin{array}{cc}
#   x_1 & = & \frac{1}{a_{11}} \left(b_1 - a_{12}x_2 - a_{13}x_3 - \cdots - a_{1n}x_n \right)\\
#   x_2 & = & \frac{1}{a_{22}} \left(b_2 - a_{21}x_1 - a_{23}x_3 - \cdots - a_{2n}x_n \right)\\
#     & \vdots & \\
#   x_n & = & \frac{1}{a_{nn}} \left(b_n - a_{n1}x_1 - a_{n2}x_2 - \cdots - a_{n,n-1}x_{n,n-1}\right)\\
# \end{array}
# $$
# 
# 2. Buat tebakan awal untuk solusinya, dinotasikan dengan $x_1^{(0)}, x_2^{(0)}, \cdots, x_n^{(0)}$ dan subtitusikan ke dalam persamaan pada langkah 1 untuk mendapatkan solusi pada iterasi pertama.
# 
# 3. Lakukan proses langkah 2 secara terus menerus sampai mendapatkan barisan solusi $x_1^{(k)}, x_2^{(k)}, \cdots, x_n^{(k)}$, untuk $k=1, 2, \cdots$, yang konvergen ke solusi dari SPL tersebut. 
# 
# 

# In[1]:


import numpy as np
from scipy.linalg import norm
import matplotlib.pyplot as plt


# __Metode Jacobi__

# In[2]:


def jacobi(A, b, x0, epsilon=1e-5, N=1000):
  k = 1
  x = np.zeros(len(b))
  while (k < N):
    for i in range(len(b)):
      sum = 0
      for j in range(len(b)):
        if i != j:
          sum += A[i,j] * x0[j]
      
      x[i] = 1/A[i,i]*(-sum + b[i])

    if norm(x - x0) < epsilon:
      break
    
    x0 = x
    k = k + 1
  return k,x


# In[3]:


def jacobi1(A, b, x0, epsilon=1e-5, N=1000):
  x = np.zeros(len(b))
  T = A - np.diag(np.diagonal(A))

  for k in range(N):
    x = (b - np.dot(T,x0))/np.diagonal(A)

    if norm(np.dot(A,x) - b) < epsilon:
      break
      
    x0 = x

  return k,x


# __Metode Gauss-Seidel__

# In[4]:


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


# __Metode SOR__

# In[5]:


def SOR(A, b, x0, omega, epsilon=1e-5, N=1000):
  x = np.zeros_like(b, dtype=np.double)
  
  for k in range(N):

    for i in range(len(b)):
      U = np.dot(A[i,:i], x[:i])
      V = np.dot(A[i,(i+1):], x0[(i+1):])
      x[i] = 1/A[i,i] * (b[i] - U - V)
      x[i] = np.dot(x0[i], (1-omega)) + np.dot(x[i], omega)
      #print(k, x) 
      
    if norm(np.dot(A,x) - b) < epsilon:
      break
    #if np.linalg.norm((x - x0)) <= epsilon:
      #break

    x0 = x

  return k,x


# __Contoh 1__:
# 

# In[6]:


A = np.array([[3, -1, 1], [3, 6, 2], [3, 3, 7]])
b = np.array([1, 0, 4])
x0 = np.array([1., 1., 1.])

x_J = jacobi1(A, b, x0)
x_GS = gauss_seidel(A, b, x0)
x_SOR = SOR(A, b, x0, 0.3)
print(x_J)
print(x_GS)
print(x_SOR)


# __Contoh 2__:
# 

# In[7]:


A = np.array([[4, 3, 0], [3, 4, -1], [0, -1, 4]])
b = np.array([24, 30, -24])
x0 = np.array([1., 1., 1.])

x_J = jacobi1(A, b, x0)
x_GS = gauss_seidel(A, b, x0)
x_SOR = SOR(A, b, x0, 0.2)
print(x_J)
print(x_GS)
print(x_SOR)


# __Contoh 4__:
# 
# Diberikan matriks
# 
# $$
# A = 
# \begin{bmatrix}
#   2 & -1 & 0 \\
#   -1 & 2 & -1 \\
#   0 & -1 & 2 \\
# \end{bmatrix}
# $$
# 
# Ingin: Melihat nilai norm matriks dari bentuk iteratif matriks $T_j, T_{gs}, T_{\omega}$
# 
# $$A = D - L - U$$
# 
# $$
# A = 
# \begin{bmatrix}
#   2 & -1 & 0 \\
#   -1 & 2 & -1 \\
#   0 & -1 & 2 \\
# \end{bmatrix}
# =
# \begin{bmatrix}
#   2 & 0 & 0 \\
#   0 & 2 & 0 \\
#   0 & 0 & 2 \\
# \end{bmatrix}
# -
# \begin{bmatrix}
#   0 & 0 & 0 \\
#   1 & 0 & 0 \\
#   0 & 1 & 0 \\
# \end{bmatrix}
# -
# \begin{bmatrix}
#   0 & 1 & 0 \\
#   0 & 0 & 1 \\
#   0 & 0 & 0 \\
# \end{bmatrix}
# $$
# 
# $T_j = D^{-1} (L+U)$, $T_{gs} = (D-L)^{-1} U$, $T_{\omega} = (D-\omega L)^{-1} [(1-\omega)D + \omega U]$

# In[8]:


D = np.array([[2.,0.,0.], [0.,2.,0.], [0.,0.,2.]])
L = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.]])
U = np.array([[0.,1.,0.], [0.,0.,1.], [0.,0.,0.]])

D_invers = np.linalg.inv(D)
LU = L+U
Tj = np.matmul(D_invers,LU)
Tj


# In[9]:


DL_invers = np.linalg.inv((D-L))
Tgs = np.matmul(DL_invers,U)
Tgs


# In[10]:


w = 1.25
A = np.linalg.inv(D-w*L)
A


# In[11]:


B = (1-w)*D + w*U
B


# In[12]:


Tsor = np.matmul(A,B)
Tsor


# Norm $l_1, l_2, l_{\infty}$ dari $T_j, T_{gs}, T_{SOR}$

# In[13]:


l = np.zeros((3,3), dtype=np.double)
l[0,0] = norm(Tj,1)
l[1,0] = norm(Tgs,1)
l[2,0] = norm(Tsor,1)

l[0,1] = norm(Tj,2)
l[1,1] = norm(Tgs,2)
l[2,1] = norm(Tsor,2)

l[0,2] = norm(Tj,np.inf)
l[1,2] = norm(Tgs,np.inf)
l[2,2] = norm(Tsor,np.inf)

l


# In[159]:




