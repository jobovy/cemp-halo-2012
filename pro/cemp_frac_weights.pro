;+
;   NAME:
;      cemp_frac_weights
;   PURPOSE:
;      calculate the fraction of cemp in inner/outer halo by using the weights
;   INPUT:
;      feh - 0 for -1.5, 1 for -2
;   OUTPUT:
;      prints result
;   HISTORY:
;      2010-04-01 - Written - Bovy (NYU)
;-
FUNCTION JACKKNIFE_FRAC_WEIGHTS, data, weight, cempcut=cempcut
IF ~keyword_set(cempcut) THEN cempcut= 1.
ndata= n_elements(data)
jack_estimates= dblarr(ndata)
FOR ii=0L, ndata-1 DO BEGIN
    IF ii EQ 0 THEN BEGIN
        thisdata= data[1:ndata-1]
        thisweight= weight[1:ndata-1]
    ENDIF ELSE IF ii EQ ndata-1 THEN BEGIN
        thisdata= data[0:ndata-2]
        thisweight= weight[0:ndata-2]
    ENDIF ELSE BEGIN
        thisdata= [data[0:ii-1],data[ii+1:ndata-1]]
        thisweight= [weight[0:ii-1],weight[ii+1:ndata-1]]
    ENDELSE
    jack_estimates[ii]= total(thisweight[where(thisdata GT cempcut)])/total(thisweight)*100.
ENDFOR
RETURN, stddev(jack_estimates)*(ndata-1.)/sqrt(double(ndata))
END
PRO CEMP_FRAC_WEIGHTS, feh=feh
IF ~keyword_set(feh) THEN feh= 0
IF feh EQ 0 THEN BEGIN
    READ_CFE_INNEROUTER, vphi, cfe, p_in, p_out, datafilename='../data/Test_for_Jo_-1.5.txt',nhdr=2,ndata=975
ENDIF ELSE BEGIN
    READ_CFE_INNEROUTER, vphi, cfe, p_in, p_out, datafilename='../data/Test_for_Jo.txt',nhdr=1,ndata=127
ENDELSE

;;Fraction of cemp in inner
frac_in= total(p_in[where(cfe GT 1.)])/total(p_in)
frac_out= total(p_out[where(cfe GT 1.)])/total(p_out)

;;Jackknife errors
efrac_in= jackknife_frac_weights(cfe,p_in)
efrac_out= jackknife_frac_weights(cfe,p_out)

print, "Inner halo: ", frac_in*100., efrac_in
print, "Outer halo: ", frac_out*100., efrac_out

END
