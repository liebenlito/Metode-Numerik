#!/usr/bin/env python
# coding: utf-8

# # Faktorisasi Matriks
# 
# Pada prinsipnya Eliminasi Gauss adalah teknik yang secara umum digunakan untuk menyelesaikan Sistem Persamaan Linear (SPL). Namun, waktu komputasi yang dibutuhkan oleh Eliminasi Gauss cenderung lama seiring bertambahnya ukuran matriks atau dapat dinotasikan $O(n^3 /3)$ {cite}`burden_numerical_2010`. Hal ini banyak dipengaruhi oleh operasi aritmatika dari proses eliminasi yang terjadi di dalam Eliminasi Gauss. Dengan demikian, kita memerlukan suatu cara agar dapat meminimalkan waktu komputasi, salah satunya dengan melakukan __faktorisasi matriks__.

# ## Faktorisasi LU
# 
# Diberikan suatu SPL yang dapat dibentuk dalam bentuk matriks-vektor sebagai berikut
# 
# $$
# A\vec{x} = \vec{b}
# $$
# 
# dimana $A \in \mathbb{R}^{n \times n}$ dan $\vec{x}, \vec{b} \in \mathbb{R}^n$. Kemudian untuk menyelesaikan SPL tersebut, kita gunakan eliminasi Gauss tanpa pivoting dan kita asumsikan juga bahwa elemen pivot di $A$ tak nol seiring berjalannya proses eliminasi pada setiap iterasinya, atau dapat ditulis $a_{ii}^{(k)}$, untuk $k = 1,2, ..., n-1$ dan $i = 1,2,...,n$. Proses eliminasi pada iterasi pertama dilakukan dengan cara
# 
# $$
# (E_j - m_{j1}E_1) \rightarrow (E_j), \hspace{1.5em} m_{j1} = \frac{a_{j1}^{(1)}}{a_{11}^{(1)}}
# $$
# 
# atau proses ini dapat direpresentasikan sebagai matriks 
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
# Proses eliminasi ini dapat juga dipandang sebagai operasi perkalian matriks antara

# In[ ]:




