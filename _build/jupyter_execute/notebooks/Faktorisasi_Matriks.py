#!/usr/bin/env python
# coding: utf-8

# # Dekomposisi Matriks
# 
# Pada prinsipnya Eliminasi Gauss adalah teknik yang secara umum digunakan untuk menyelesaikan Sistem Persamaan Linear (SPL). Namun, waktu komputasi yang dibutuhkan oleh Eliminasi Gauss cenderung lama seiring bertambahnya ukuran matriks atau dapat dinotasikan $O(n^3 /3)$ {cite}`burden_numerical_2010`. Hal ini banyak dipengaruhi oleh operasi aritmatika dari proses eliminasi yang terjadi di dalam Eliminasi Gauss. Dengan demikian, kita memerlukan suatu cara agar dapat meminimalkan waktu komputasi, salah satunya dengan melakukan __dekomposisi matriks__.

# ## Dekomposisi LU
# 
# Diberikan suatu SPL yang dapat dibentuk dalam bentuk matriks-vektor sebagai berikut
# 
# $$
# A\vec{x} = \vec{b}
# $$
# 
# dimana $A \in \mathbb{R}^{n \times n}$ dan $\vec{x}, \vec{b} \in \mathbb{R}^n$. Kemudian untuk menyelesaikan SPL tersebut, kita gunakan eliminasi Gauss tanpa pivoting dan kita asumsikan juga bahwa elemen pivot di $A$ tak nol seiring berjalannya proses eliminasi pada setiap iterasinya, atau dapat ditulis $a_{ii}^{(k)}$, untuk $k = 1,2, ..., n$ dan $i = 1,2,...,n$. Proses eliminasi pada iterasi pertama dilakukan dengan cara
# 
# $$
# (E_j - m_{j1}E_1) \rightarrow (E_j), \hspace{1.5em} m_{j1} = \frac{a_{j1}^{(1)}}{a_{11}^{(1)}}
# $$
# 
# atau operasi ini dapat direpresentasikan sebagai __matriks transformasi Gaussian pertama__
# 
# $$
# M^{(1)} = 
# \begin{bmatrix}
# 1       & 0 & \cdots & 0 \\
# -m_{21} & 1 &        & 0 \\
# \vdots  &   & \ddots & \vdots \\
# -m_{n1} & 0 & \cdots & 1 \\
# \end{bmatrix}
# $$
# 
# Sehingga, hasil perkalian matriks antara $M^{(1)}$ dan $A^{(1)}$ adalah
# 
# $$
# A^{(2)}\vec{x} = M^{(1)}A\vec{x} = M^{(1)} \vec{b} = \vec{b}^{(2)}.
# $$
# 
# Secara umum, ketika kita melakukan perkalian matriks A dengan matriks transformasi Gaussian pada iterasi ke $k$ akan menghasilkan
# 
# $$
# A^{(k+1)} \vec{x} = M^{(k)} A^{(k)} \vec{x} = M^{(k)} \cdots M^{(1)} A \vec{x} = M^{(k)}\vec{b}^{(k)} = \vec{b}^{(k+1)} = M^{(k)} \cdots M^{(1)} \vec{b}.
# $$
# 
# Proses ini akan berhenti pada $A^{(n)} \vec{x} = \vec{b}^{(n)}$, dimana $A^{(n)}$ menjadi matriks segitiga atas
# 
# $$
# A^{(n)} = 
# \begin{bmatrix}
# a_{11}^{(1)} & a_{12}^{(1)} & \cdots & a_{1n}^{(1)} \\
# 0            & a_{22}^{(2)} &        & a_{2n}^{(2)} \\
# \vdots       &   & \ddots & \vdots \\
# 0            & 0 & \cdots & a_{nn}^{(n)} \\
# \end{bmatrix}
# $$
# 
# diberikan oleh 
# 
# $$
# A^{(n)} = M^{(n-1)} M^{(n-1)} \cdots M^{(1)} A.
# $$

# Kemudian kita panggil kembali 
# 
# $$
# A^{(k+1)} \vec{x} = M^{(k)} A^{(k)} \vec{x} = M^{(k)} \cdots M^{(1)} A \vec{x} = M^{(k)}\vec{b}^{(k)} = \vec{b}^{(k+1)}
# $$
# 
# dimana $M^{(k)}$ dihasilkan dari operasi baris
# 
# $$
# (E_j - m_{jk}E_k) \rightarrow (E_j), \hspace{1.5em} j=k+1, \cdots, n.
# $$
# 
# Jika kita ingin mengembalikan nilai dari matriks $A^{(k+1)} ke A^{(k)}$, maka kita perlu operasi balikan 
# 
# $$
# (E_j + m_{jk}E_k) \rightarrow (E_j), \hspace{1.5em} j=k+1, \cdots, n,
# $$
# 
# berarti ini sama saja mencari invers dari matriks $M$, $[M^{(k)}]^{-1}$. Dengan demikian, jika kita teruskan operasi invers ini sampai dengan iterasi pertama, maka akan menghasilkan matriks $A$, dan operasi invers dari matriks $M$ dapat kita tuliskan ekivalen dengan matriks segitiga bawah $L$ atau dapat ditulis
# 
# $$
# L = L^{(1)} L^{(2)} \cdots L^{(n-1)} = 
# \begin{bmatrix}
# 1       & 0 & \cdots & 0 \\
# m_{21} & 1 &        & 0 \\
# \vdots  &   & \ddots & \vdots \\
# m_{n1} & \cdots & m_{n,n-1} & 1 \\
# \end{bmatrix}
# $$
# 
# Jadi, $A$ dapat kita dekomposisi menjadi hasil kali dari matriks segitiga bawah $L$ dan matriks segitiga atas $U$, atau dapat ditulis
# 
# $$
# \begin{array}{lll}
# LU &=& L^{(1)}L^{(2)} \cdots L^{(n-3)} L^{(n-2)} L^{(n-1)} . M^{(n-1)} M^{(n-2)} M^{(n-3)} \cdots M^{(2)} M^{(1)} A \\
# &=& [M^{(1)}]^{-1} [M^{(2)}]^{-1} \cdots [M^{(n-2)}]^{-1} [M^{(n-1)}]^{-1} . M^{(n-1)} M^{(n-2)} M^{(n-3)} \cdots M^{(2)} M^{(1)} A \\
# &=& A.
# \end{array}
# $$

# ```{admonition} Teorema
# Jika Eliminasi Gauss tanpa penukaran baris dapat digunakan untuk menyelesaikan $A\vec{x} = \vec{b}$, maka matriks $A$ dapat didekomposisi ke dalam perkalian matriks antara matriks segitiga bawah $L$ dan matriks segitiga atas $U$, $A = LU$
# 
# $$
# U = 
# \begin{bmatrix}
# a_{11}^{(1)} & a_{12}^{(1)} & \cdots & a_{1n}^{(1)} \\
# 0            & a_{22}^{(2)} &        & a_{2n}^{(2)} \\
# \vdots       &   & \ddots & \vdots \\
# 0            & 0 & \cdots & a_{nn}^{(n)} \\
# \end{bmatrix}
# , L = 
# \begin{bmatrix}
# 1       & 0 & \cdots & 0 \\
# m_{21} & 1 &        & 0 \\
# \vdots  &   & \ddots & \vdots \\
# m_{n1} & \cdots & m_{n,n-1} & 1 \\
# \end{bmatrix}
# ,
# $$
# 
# dimana $m_{ji} = \frac{a_{ji}^{(i)}}{a_{ii}^{(i)}}$, untuk $i,j = 1,2,\cdots, n$.
# ```

# In[1]:


import numpy as np


# In[2]:


def faktorisasiLU(A):
    n = len(A)
    L = np.eye(n, dtype=np.double)

    for i in range(n):
        faktor = A[i+1:, i]/A[i,i]
        L[i+1:, i] = faktor
        A[i+1:] -= faktor[:, np.newaxis]*A[i]
    
    U = A
    
    return L,U


# In[3]:


A = np.array([[1, 1, 0, 3],
              [2, 1, -1, 1],
              [3, -1, -1, 2],
              [-1, 2, 3, -1]], dtype=np.float64)

b = np.array([4, 1, -3, 4], dtype=np.float64)


# In[4]:


def subtitusi_maju(L,b):
    n = len(b)
    L_c = np.c_[L,b]
    y = np.zeros_like(b, dtype=np.double)

    # Lakukan subtitutsi mundur
    y[0] = b[0]/L[0,0]
    for i in range(1,n):
        y[i] = (b[i] - np.dot(L[i,:i], y[:i])) / L[i,i]
    
    return y


# In[5]:


def subtitusi_mundur(U,y):
    n = len(y)
    x = np.zeros_like(y, dtype=np.double)
    
    # Lakukan subtitutsi mundur
    x[-1] = y[-1]/U[-1,-1]
    for i in range(n-2, -1, -1):
        x[i] = (y[i] - np.dot(U[i,i:], x[i:])) / U[i,i]

    return x


# In[6]:


def lu(A,b):
    L, U = faktorisasiLU(A)
    y = subtitusi_maju(L,b)

    return subtitusi_mundur(U,y)


# In[7]:


x = lu(A,b)
x


# In[ ]:




