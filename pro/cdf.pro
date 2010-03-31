;+
;   NAME:
;      cdf
;   PURPOSE:
;      calculate the CDF of a data set
;   CALLING SEQUENCE:
;      cdf, data, xvec, cdf, weight=weight, /invert
;   INPUT:
;      data - the data (array)
;   OPTIONAL INPUT:
;      weight - weights to be applied to the data
;   KEYWORDS:
;      invert - calculate 1-CDF
;   OUTPUT:
;      xvec - xs in which the cdf is calculated
;      cdf - the cdf (can be plotted as plot, xvec, cdf
;   HISTORY:
;      2010-03-31 - Written - Bovy (NYU)
;-
PRO CDF, data, xvec, cdf, weight=weight,invert=invert
ndata= n_elements(data)
IF ~keyword_set(weight) THEN weight= dblarr(ndata)+1.
xvec= dblarr(ndata)
cdf= dblarr(ndata)
IF keyword_set(invert) THEN sortindx= REVERSE(SORT(data)) ELSE sortindx= SORT(data)
FOR ii=0L, ndata-1 DO BEGIN
    xvec[ii]= data[sortindx[ii]]
    cdf[ii]= total(weight[sortindx[0:ii]])
ENDFOR
cdf= cdf/cdf[ndata-1]
END
