Title
=====

# sub-title

+++ this is a block break

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

