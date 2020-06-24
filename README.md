# Welcome
This is Fortran-90 implementation of ECGMM algorithm presented in:
Hao, J., Koester, B. P., Mckay, T. a., Rykoff, E. S., Rozo, E., Evrard, A., Annis, J., et al. (2009). Precision Measurements of the Cluster Red Sequence Using an Error-Corrected Gaussian Mixture Model. The 
Astrophysical Journal, 702(1), 745–758. [ADS record](http://adsabs.harvard.edu/abs/2009ApJ...702..745H)

Go ahead and try:

```
$ git clone https://github.com/minzastro/ecgmm_fortran
```

It is designed to be used with Python via f2py.

Here's an example of some Python code:

```
#!python
#!/usr/bin/python
import numpy as np
from ecgmm import ecgmm
from numpy.random import normal

data = np.hstack((normal(2.0, 0.1, 150), normal(1.5, 3.0, 30)))
err = np.zeros(len(data))
w = np.array([0.75, 0.25])
mu = np.array([2., 0.0])
sigma = np.array([1.0, 2.0])
ecgmm.iterate(data, err, w, mu, sigma, 50)
print w[0], w[1], mu[0], mu[1], sigma[0], sigma[1]
```
