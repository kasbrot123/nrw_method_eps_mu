# Nicholson-Ross-Weir Method (Python)

Nicholson-Ross-Weir (NRW) method in python for calculating refractive index,
permittivity and permeability from scattering parameters (S11/S21). The
algorithm is also taking into account the phase ambiguity of the transmission
coefficient and selects the best n based on the group velocity.

## Content

- **Function for NRW Method**: Calculates the permittivity and permeability.
  Takes frequency values, complex S11/S21, the sample length and a cutoff
  frequency (*optional*).
- **Function for reverse NRW Method**: Calculates the scattering parameters
  S11/S21 from the permittivity and the permeability. Takes frequency values,
  complex permittivity/permeability, the sample length and a cutoff frequency
  (*optional*).
- **High precision forwards-backwards NRW calculation**: This was done as an
  instability analysis. One can specify the permittivity/permeability to
  calculate the scattering parameters S11/S21 and back to
  permittivity/permeability using NRW method. By introducing noise it is
  possible to test a measurement result. 


## Analysis

### Basis

The NRW method starts using the complex scattering parameters $S_{11}$ and
$S_{21}$. The reflection coefficient is defined as:

```math
X = \frac{S_{11}^2-S_{21}^2+1}{2S_{11}}
```
```math
\Gamma =  X \pm \sqrt(X^2-1)
```
Choose the sign of the root for $\Gamma$ such that $|\Gamma|\le1$.
The transmission coefficient is defined as.

```math
T = \frac{S_{11}+S_{21}-\Gamma}{1-(S_{11}+S_{21})\Gamma}
```

Using $\Gamma$ and $T$ one can calculate $\Lambda$:

```math
\frac{1}{\Lambda^2} = (\frac{\varepsilon_r\mu_r}{\lambda_0}-\frac{1}{\lambda_c^2}) = -(\frac{1}{2\pi L}ln(\frac{1}{T}))^2
```
The equation for $\Lambda$ is not well defined since the logarithm of $1/T$ has
multiple solutions which are equal to $ln(1/T) + i(\Theta + 2\pi n)$ where $n$
is an integer value. In the NRW method one can estimate the value $n$ by using
the group delay. 

The idea is that the group delay is the derivative of the phase and therefore
is independent of the $2\pi n$. One can compare the measured group delay and
the calculated group delay and the best value for $n$ is found where:

```math
\tau_{meas} - \tau_{calc} \approx 0
```

The group delay for both values is defined as:

```math
\tau = \frac{d\phi}{\omega} = -\frac{1}{2\pi}\frac{d\phi}{df}
```
```math
\tau_{meas} = -\frac{1}{2\pi}\frac{d\phi_{meas}}{df} = -\frac{1}{2\pi}\frac{d}{df}arg(T)
```
```math
\tau_{calc} = -\frac{1}{2\pi}\frac{d\phi_{calc}}{df} = \frac{d}{df}\frac{L}{\Lambda} = 
```
```math
=L\frac{d}{df}\sqrt{\frac{\varepsilon_r\mu_r f^2}{c^2}-\frac{1}{\lambda_c^2}} = \frac{1}{c^2}\frac{f\varepsilon_r\mu_r+f^2\frac{1}{2}\frac{d(\varepsilon_r\mu_r)}{df}}{\sqrt{\frac{\varepsilon_r\mu_r f^2}{c^2}-\frac{1}{\lambda_c^2}}}L
```


With this the permittivity and permeability is defined:

```math
\mu_r = \frac{1+\Gamma_1}{\Lambda(1-\Gamma)\sqrt{\frac{1}{\lambda_0^2}-\frac{1}{\lambda_c^2}}}
```
```math
\varepsilon_r = \frac{\lambda_0^2}{\mu_r}\left(\frac{1}{\lambda_c^2} - \left[\frac{1}{2\pi L}\ln\left(\frac{1}{T}\right)\right]^2\right)
```


### Literature

The analysis explained in literature.

- O. Luukkonen, S. I. Maslovski and S. A. Tretyakov, "A Stepwise
  Nicolson–Ross–Weir-Based Material Parameter Extraction Method," in IEEE
  Antennas and Wireless Propagation Letters, vol. 10, pp. 1295-1298, 2011, doi:
  10.1109/LAWP.2011.2175897.
- A. M. Nicolson and G. F. Ross, "Measurement of the Intrinsic Properties of
  Materials by Time-Domain Techniques," in IEEE Transactions on Instrumentation
  and Measurement, vol. 19, no. 4, pp. 377-382, Nov. 1970, doi:
  10.1109/TIM.1970.4313932.
- W. B. Weir, "Automatic measurement of complex dielectric constant and
  permeability at microwave frequencies," in Proceedings of the IEEE, vol. 62,
  no. 1, pp. 33-36, Jan. 1974, doi: 10.1109/PROC.1974.9382.
- Measurement of Dielectric Material Properties Application Note R&S.



## Instabilities

Some instabilities were observed when:

- The cutoff frequency is near the lowest measurement frequency
- Long sample lengths or high permittivity/permeability values (short
  wavelengths media)

