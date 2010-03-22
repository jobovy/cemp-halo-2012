;+
;   NAME:
;      test_weighted_kstwo
;   PURPOSE:
;      simple testing code for weighted_kstwo: draw samples from a
;      mixture of two Gaussians, calculate weights, perform the KS
;      test with and without weights
;   INPUT: (optional)
;      seed - seed for the random number generator
;      n1 - number of samples to draw for the first data set (one
;           Gaussian)
;      n2 - number of samples to draw for the second data set (mixture
;           of two Gaussians, first one equal to that of data1)
;      plot - plot histograms to output
;   OUTPUT:
;   HISTORY:
;      2010-03-22 - Written - Bovy (NYU)
;-
PRO TEST_WEIGHTED_KSTWO, seed=seed, n1=n1, n2=n2, plot=plot
IF ~keyword_set(seed) THEN seed= -1L
IF ~keyword_set(n1) THEN n1= 100
IF ~keyword_set(n2) THEN n2= 100
data1= randomn(seed,n1)
mean2= 5.
sig2= 2.
data2= dblarr(n2)
weight2= dblarr(n2)
FOR ii=0L, n2-1 DO BEGIN
    onequestionmark= randomu(seed)
    if onequestionmark LT 0.5 then data2[ii]= randomn(seed) else data2[ii]= randomn(seed)*sig2+mean2
    weight2[ii]= exp(-0.5*data2[ii]^2.)/(exp(-0.5*data2[ii]^2.)+1./sig2*exp(-0.5*(data2[ii]-mean2)^2./sig2^2.))
ENDFOR

print, "With weights:"
weighted_kstwo, data1, data2, D, prob, weight2=weight2 & print, prob
weight1= dblarr(n1)+1.0/n1
print, "Without weights:"
weighted_kstwo, data1, data2, D, prob & print, prob
weight= dblarr(n2)+1./n2
IF keyword_set(plot) THEN BEGIN
    print, "White: data2 with weights, Red: data1 (compare red and white), blue: data2 without weights"
    hogg_plothist, data2, weight=weight2/total(weight2), /totalweight, xvec=xvec
    npix=n_elements(xvec)
    hogg_plothist, data1, weight=weight1, /totalweight, /overplot, color=djs_icolor('red'), npix=npix
    hogg_plothist, data2, /overplot, weight=weight, color=djs_icolor('blue'), /totalweight, npix=npix
ENDIF
END
