This code implements the measurement of component overlapping in the Gaussian mixture model, \textbf{Overlap Rate (OLR)}, proposed by \textit{Haojun Sun}. ([Measuring the component overlapping in the Gaussian mixture model](https://link.springer.com/article/10.1007/s10618-011-0212-3))

The OLR is defined as:
$$OLR(G_1, G_2) = 
\begin{cases}
	1 &\mbox{if $p(X)$ has one peak} \\
	\frac{p(X_{saddle})}{p(X_{submax})} &\mbox{if $p(X)$ has two peaks}
\end{cases}
$$

where $X_{saddle}$ is the saddle point of pdf and $p(X_{submax})$ is the lower peak point of pdf.

This paper proves that the saddle point and means are on the same ridge curve, thus we can find the $X_{saddle}$ along the ridge curve. The ridge curve can be described by a specific equation.

This code implements the algorithm presented in this paper in Python.

