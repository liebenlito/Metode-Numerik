Sistem Bilangan Floating-Point
=====
Bilangan di sistem floating-point direpresentasikan sebagai deret dari digit dimana digit tersebut merepresentasikan suatu angka yang berbeda. 
```{note}
Secara umum, sistem floating-point ternormalisasi dapat ditulis sebagai 

$$F = \pm d_1 . d_2 d_3 d_4 \ldots d_p \times \beta^E$$

dimana
1. $\pm$ adalah bit yang merepresentasikan tanda dari bilangan
2. $d_1 . d_2 d_3 d_4 \ldots d_t$ disebut *mantissa*.  Digit-digit $d_2 d_3 d_4 \ldots d_p$ disebut *fraction/fraksi* dengan $t$ digit presisi.  Sistem ternormalisasi secara spesifik $d_1 \neq 0$ kecuali bilangan tersebut adalah $0$.
3. $\beta$ adalah *basis*.  Untuk binary $\beta = 2$, untuk desimal $\beta = 10$, dst.
4. $E$ adalah *eksponen*, suatu bilangan bulat antara $[E_{\min}, E_{\max}]$.
```
Bagian penting dalam sistem floating-point yaitu
1. Terdapat himpunan diskrit dan berhingga yang merepresentasikan suatu bilangan
2. Dapat merepresentasikan bilangan yang tidak terdisrtibusi secara teratur pada garis bilangan (real)
3. Aritmatika di dalam sistem floating-point menghasilkan hasil yang berbeda dibandingkan dengan aritmatika di dalam sistem bilangan real

# Sifat-Sifat dari sistem floating-point
Semua sistem floating-point juga berupa beberapa bilangan yang penting:
 - Bilangan ternormalisasi terkecil (underflow)
 - Bilangan ternormalisasi terbesar (overflow)
 - Nol
 - Machine $\epsilon$ atau $\epsilon_{\text{machine}}$
 - `inf` dan `nan`, tak hingga dan **N**ot **a** **N**umber

Contoh:
Misalkan terdapat sistem desimal 2-digit (ternormalisasi)

$$f = \pm d_1 . d_2 \times 10^E$$

dengan $E \in [-2, 0]$.

**Bilangan dan distribusi dari bilangan**

1. Berapa banyak bilangan yang dapat direpresentasikan dengan sistem ini?

$$
    f = \pm d_1 . d_2 \times 10^E ~~~  E \in [-2, 0]
$$

$$ 
    2 \times 9 \times 10 \times 3 + 1 = 541
$$

2. Bagaimana distribusinya pada garis bilangan "real"?
```python
import matplotlib.pyplot as plt

nilai_d1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
nilai_d2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
nilai_E = [0, -1, -2]

fig = plt.figure(figsize=(10.0, 1.0))
ax = fig.add_subplot(1, 1, 1)

for E in nilai_E:
    for d1 in nilai_d1:
        for d2 in nilai_d2:
            ax.plot( (d1 + d2 * 0.1) * 10**E, 0.0, 'r+', markersize=20)
            ax.plot(-(d1 + d2 * 0.1) * 10**E, 0.0, 'r+', markersize=20)
            
ax.plot(0.0, 0.0, '+', markersize=20)
ax.plot([-10.0, 10.0], [0.0, 0.0], 'k')

ax.set_title("Distribusi Nilai")
ax.set_yticks([])
ax.set_xlabel("x")
ax.set_ylabel("")
ax.set_xlim([-2, 2])
plt.show()

```

3. Apa nilai dari underflow dan overflow?

Bilangan terkecil dapat direpresentasikan dengan underflow:  $1.0 \times 10^{-2} = 0.5$

Bilangan terbesar dapat direpresentasikan dengan overflow:  $1.9 \times 10^0 = 1.9$

Catatan bahwa semua sistem floating-point IEEE 754 menggunakan bilangan binari.  

Cara cepat:
$$
    2^3 2^2 2^1 2^0 . 2^{-1} 2^{-2} 2^{-3}
$$
melambangkan 8, 4, 2, 1 . 1/2, 1/4, 1/8, ...

[def]

[def]: http://www.abc.com "abc"

```{note}
Here is an *admonition* directive.
You can write **any** Markdown in here, and it will be `syntax highlighted`
```

```{tip}
Even nested directives will work:
```

```python
def f():
    pass
```

$$
\begin{align}
y &= \lim_{n \rightarrow \infty}\sum_{i=0}^{n} \frac{1}{n} \\
 &= 1
\end{align}
$$

## Review Aljabar Linear

```{admonition} Definisi 
Misalkan $\lbrace \textbf{v}^{(1)}, \textbf{v}^{(2)}, \cdots, \textbf{v}^{(k)}\rbrace$ adalah himpunan vektor. Himpunan tersebut dikatakan __linearly independent__ jika

$$
\textbf{0} = \alpha_1 \textbf{v}^{(1)} + \alpha_2 \textbf{v}^{(2)} + \cdots + \alpha_k \textbf{v}^{(k)}
$$

maka $\alpha_i = 0$ untuk setiap $i = 0, 1, \cdots, k$. Begitu pula sebaliknya, maka himpunan vektornya disebut __linearly dependent__. 
```

```{admonition} Teorema 
Diberikan $\lbrace \textbf{v}^{(1)}, \textbf{v}^{(2)}, \cdots, \textbf{v}^{(k)}\rbrace$ adalah himpunan vektor yang linearly independent di $\mathbb{R}^n$. Maka untuk setiap vektor $\textbf{x} \in \mathbb{R}^n$ 

$$
\textbf{0} = \alpha_1 \textbf{v}^{(1)} + \alpha_2 \textbf{v}^{(2)} + \cdots + \alpha_k \textbf{v}^{(k)}
$$

maka $\alpha_i = 0$ untuk setiap $i = 0, 1, \cdots, k$. Begitu pula sebaliknya, maka himpunan vektornya disebut __linearly dependent__. 
```

