# Details on predicting photometry

It is sometimes not evident that there are essential differences between photometric systems. But even less known, the difference between detector count types (energy or photons) requires special care.

This section reviews the essential details for computing the luminosity
and magnitude of a star through a photometric passband filter. We do not
discuss calibration, which instrument documentations cover, in principle.

Let's consider a filter throughput (a.k.a, transmission curve, or response function) defined in wavelength by the dimensionless function $T(\lambda)$. This function tells you what fraction of the arriving photons at wavelength $\lambda$ get through the instrument.  Therefore, the total number of
photons, per unit time per unit area, received in this filter is

$$N_{tot} = \frac{1}{hc} \int_\lambda f_\lambda\,\lambda\,T(\lambda)\,d\lambda,$$

where $f_\lambda$ is the wavelength-dependent flux density of an object
given in energy per unit time per unit area per unit wavelength.

Consequently, interpreting $\lambda T(\lambda)$ as a distribution leads to the statistical mean of the flux density, $\overline{f_\lambda}$:

$$\overline{f_\lambda}(T) = \frac{\int_\lambda \lambda f_\lambda T(\lambda) d\lambda}{\int_\lambda \lambda T(\lambda) d\lambda}.$$

Note that this is not the mean flux density because of the $\lambda$ factor in the integrals. It is the mean photon rate density in this filter commonly expressed in stellar physics literature as
$erg.s^{-1}.cm^{-2}.Å^{-1}$ or $W.m^{-2}.nm^{-1}$.

⚠️ from pyphot example:
```python
# computing the flux of a spectrum
flux = lib['hst_wfc3_f110w'].get_flux(lamb, spec)
# lamb may have units, otherwise assuming consistent definitions.

# computing the flux of many spectra
fluxes = lib['hst_wfc3_f110w'].get_flux(lamb, spectra, axis=1)
```

Finally, at least for instruments using CCD or CCD-like cameras, i.e., counting photons, we obtain the usual definition of a magnitude

$$mag_\lambda(T) = -2.5\,\log_{10}\left(\overline{f_\lambda}\right) - ZP\left(\overline{f_\lambda}\right),$$

where $ZP(\overline{f_\lambda})$ gives the pass-band reference value
(zeropoint) for a given photometric/magnitude system.

However, the zeropoints themselves depend on the adopted photometric system used to report the measurements. They may vary fundamentally from one to another.  Below we briefly describe the primary systems used in large surveys.


## Vega magnitude system

This system is defined such that the star Alpha Lyr (Vega) has magnitude 0 in any passband filter. In other words, the zeropoints are set by the magnitude of vega, $-2.5 \log_{10} \overline{f_\lambda}(Vega)$, or

$$mag_{Vega}(T) = -2.5\,\log_{10}\left(\overline{f_\lambda} / \overline{f_\lambda}(Vega)\right).$$

⚠️ from pyphot example:
```python
# convert to magnitudes
import numpy as np
f = lib['hst_wfc3_f110w']
fluxes = f.get_flux(lamb, spectra, axis=1)
mags = -2.5 * np.log10(fluxes) - f.Vega_zero_mag
# or similarly
mags = -2.5 * np.log10(fluxes / f.Vega_zero_flux)
```

## Johnson system

The Johson system is defined such that the star Alpha Lyr (Vega) has $V=0.03$ and all colors equal to zero. It is very similar to the Vega magnitude system, but **using mean flux definition** (instead of photon counts), i.e., **energy
counter** detectors

$$\widetilde{f_\lambda}(T) = \frac{\int_\lambda f_\lambda T(\lambda) d\lambda}{\int_\lambda T(\lambda) d\lambda},$$

(the true definition of mean flux throughout a given transmission filter.)

> **Note**: Table A2 of Bessell et al. 1998 gives zero points for the UBVRIJHKL(+Kp and L’) filters in the Counsins-Glass-Johnson system.

If one defines the **effective wavelength** $\lambda_{\rm eff}$ as the
photon weighted mean wavelength:

$$\lambda_{\rm eff} = \frac{\int \lambda f_\lambda T(\lambda) d\lambda}{\int f_\lambda T(\lambda) d\lambda},$$

```python
# the effective wavelength for vega is given by
lib['ground_johnson_u'].leff
```

then the difference between the Johnson and Vega systems within the same filter is given by

$$\widetilde{mag}_\lambda - \overline{mag}_\lambda = 0.03 - 2.5 \log_{10} \frac{\lambda_{\rm eff}(Vega)}{\lambda_{\rm eff}(star)}, $$
where we explicit which equation was used to compute magnitudes.


```python
# The switch between the energy and the photon count equation is done
# through the `Filter.set_dtype` method, and becomes transparent for any
# use. So if you define you own filter either use the constructor or the
# method

# define a constant filter in energy count from 100 to 110 AA
f = Filter(np.arange(100, 110), np.ones(10), \
           dtype='energy', unit='AA')
# manually set the detector type
f.set_dtype('photon')
```



## AB magnitude system

This system is defined such that, when monochromatic flux $f_\nu$ is measured in $erg\,s^{-1}\,cm^{-2} Hz^{-1}$,

$$mag_{AB}(T) = -2.5\, \log_{10}(\overline{f_\nu}) - 48.60,$$
where the value of the constant is selected to define $m_{AB}=V$ for a
flat-spectrum source. In this system, an object with constant flux per unit frequency interval has zero color.

Koornneef et al. gives the respective definition of $\overline{f_\nu}(T)$:

$$\overline{f_\nu}(T) = \frac{\int_\nu f_\nu T(\nu) d\nu / \nu}{\int_\nu T(\nu) d\nu / \nu}
= \frac{\int_\lambda f_\nu T(\lambda) d\lambda / \lambda}{\int_\lambda T(\lambda) d\lambda / \lambda}$$

To go back to wavelength units, we have $d\nu = (c/\lambda^2) d\lambda$.

If one defines the **pivot wavelength** $\lambda_p$ to convert between
$\overline{f_\nu}$ and $\overline{f_\lambda}$ as

$$\overline{f_\nu} = \frac{\lambda_p^2}{c} \overline{f_\lambda},$$

one can easily show that

$$\lambda_p^2 = \frac{\int_\lambda T(\lambda)\,\lambda\,d\lambda}{\int_\lambda T(\lambda)\,d\lambda /\lambda}.$$

Therefore for filters with AB magnitudes, one can compute

$$mag_{AB}(T) = -2.5\, \log_{10}(\overline{f_\lambda}) - 2.5\log_{10}\left(\lambda_p^2/c\right) - 48.6,$$
where one must care to use the speed of light $c$ and $\lambda_p$ in matching units.

```python
# convert to magnitudes
import numpy as np
f = lib['hst_wfc3_f110w']
fluxes = f.get_flux(lamb, spectra, axis=1)
mags = -2.5 * np.log10(fluxes) - f.AB_zero_mag
# or similarly
mags = -2.5 * np.log10(fluxes / f.AB_zero_flux)
```



## ST magnitude system

This system is defined such as a source with flat $f_\lambda$ will have the same magnitude in every filter.

Koornneef et al. (1986; same as above) defines
$$mag_{ST}(T) = -2.5\, \log_{10}(\overline{f_\lambda}) - 21.1,$$

```python
# convert to magnitudes
import numpy as np
f = lib['hst_wfc3_f110w']
fluxes = f.get_flux(lamb, spectra, axis=1)
mags = -2.5 * np.log10(fluxes) - f.ST_zero_mag
# or similarly
mags = -2.5 * np.log10(fluxes / f.ST_zero_flux)
```


## Jansky definition

The jansky (symbol $Jy$) is a non-SI unit of spectral flux density, it is equivalent to $10^{−26} W.m^{-2}.Hz^{-1}$ or $10^{-23} erg.s^{-1}.cm^{-2}.Å^{-1}$.

```python
import numpy as np
f = lib['hst_wfc3_f110w']
print(f.AB_zero_Jy, f.Vega_zero_Jy, f.ST_zero_Jy)
```

---
## References

* Bessel, M. S. 1990, PASP, 91, 589;

* Bessel, M. S. 1983, PASP, 95, 480;

* Bessel, M. S. 1990, PASP, 102, 1181;

* Hayes, D. S., \& Latham, D. W. 1975, ApJ, 197, 593;

* Johnson, H. L. \& Morgan, W. W. 1953, ApJ, 117, 313

* Oke, J.B. 1974, ApJS, 27, 21;

* Koornneef, Bohlin, Buser, Horne, Turnshek: Synthetic photometry and the calibration of HST.
