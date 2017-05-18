# GMM-overlap-rate

This code implements the measurement of component overlapping in the Gaussian mixture model, **Overlap Rate (OLR)**, proposed by *Haojun Sun*. ([Measuring the component overlapping in the Gaussian mixture model](https://link.springer.com/article/10.1007/s10618-011-0212-3))

The OLR is the ratio to the saddle point of pdf of a Gaussian mixture distribution to the lower peak point of pdf of this distribution.

This paper proves that the saddle point and means are on the same ridge curve, thus we can find the $X_{saddle}$ along the ridge curve. The ridge curve can be described by a specific equation.

This code implements the algorithm presented in this paper in Python.


