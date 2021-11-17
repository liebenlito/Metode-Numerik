#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/liebenlito/Metode-Numerik/blob/main/model_SIR.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# # Model SIR 
# Misalkan di dalam suatu pengamatan pada waktu ke $t$ di dalam populasi __tertutup__, terdapat subpopulasi yang rentan (Susceptible, $S$), subpopulasi yang terjangkit penyakit (Infection, $I$), dan subpopulasi yang sembuh (Recovery, $R$). Model populasi ini dinamakan model $SIR$. Model SIR memiliki asumsi-asumsi sebagai berikut:
# 
# 1. Individu yang sembuh mempunyai kekebalan terhadap penyakit tersebut (tidak bisa menjadi rentan lagi).
# 2. Penyakit tidak fatal (tidak menyebabkan kematian).
# 3. Masa inkubasi sangat singkat. Masa inkubasi adalah masa antara mulai terjadinya infeksi sampai dengan tampaknya gejala.
# 4. Laju kenaikan proporsi (jumlah individu) dari kelas terjangkit dan menularkan sebanding dengan proporsi (jumlah individu) dari kelas rentan dan kelas terjangkit. Proporsi (jumlah individu) kelas rentan minimum dengan laju yang sama, yaitu $\beta SI$.
# 5. Laju kenaikan proporsi (jumlah individu) kelas sembuh sebanding dengan proporsi ( jumlah individu) dari kelas terjangkit dan menularkan, yaitu $\alpha I$.
# 
# 
# 
# 
# 
# 
# 

# In[ ]:




