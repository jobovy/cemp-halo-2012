;+
;   NAME:
;      feh_cemp
;   PURPOSE:
;      derive and fit the cemp vs. feh relation
;   INPUT:
;      hardcut - if set, use a hardcut of pij=0.7
;   OUTPUT:
;   HISTORY:
;      2010-06-01 - Written - Bovy (NYU)
;-
FUNCTION JACKKNIFE_FRAC_HARDCUT, data, cempcut=cempcut
IF ~keyword_set(cempcut) THEN cempcut= 1.
ndata= n_elements(data)
jack_estimates= dblarr(ndata)
FOR ii=0L, ndata-1 DO BEGIN
    IF ii EQ 0 THEN BEGIN
        thisdata= data[1:ndata-1]
    ENDIF ELSE IF ii EQ ndata-1 THEN BEGIN
        thisdata= data[0:ndata-2]
    ENDIF ELSE BEGIN
        thisdata= [data[0:ii-1],data[ii+1:ndata-1]]
    ENDELSE
    jack_estimates[ii]= double(n_elements(thisdata[where(thisdata GT cempcut)]))/n_elements(thisdata)*100.
ENDFOR
RETURN, stddev(jack_estimates)*(ndata-1.)/sqrt(double(ndata))
END
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
PRO FEH_CEMP, hardcut=hardcut
datafilename= '../data/File_Test_Carbon.txt'
readcol, datafilename, blah, feh, V, sV, cfe, format='A,D,D,D,D'

fehmin= -1.4
fehstep= -0.2
fehs= dblarr(7)
cemps= dblarr(7)
FOR ii=0L, 5 DO BEGIN
    fehs[ii]= fehmin+(ii+0.5)*fehstep
    cemps[ii]= n_elements(where(feh LT (fehmin+ii*fehstep) AND $
                            feh GE (fehmin+(ii+1)*fehstep) AND $
                            cfe GE 0.9))/double(n_elements(where(feh LT (fehmin+ii*fehstep) AND $
                            feh GE (fehmin+(ii+1)*fehstep))))
ENDFOR
fehs[6]= -2.8
cemps[ii]= n_elements(where(feh LT -2.6 AND $
                            feh GE -3 AND $
                            cfe GE 0.9))/double(n_elements(where(feh LT -2.6 AND $
                                                                 feh GE -3.)))

djs_plot, fehs, cemps, psym=5, xrange=[-1.4,-3.5], yrange=[0,1], $
  xtitle='[Fe/H]', ytitle='CEMP fraction'

coeffs= poly_fit(fehs,cemps,2)
xs= dindgen(1001)/1000*(-1.6)-1.4
ys= coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.
djs_oplot, xs, ys, linestyle=2

indx= where(feh LT -2. and feh GT -2.5)
amp= [0.5,0.5]
xmean= [0.,-100.]
xcovar=[1000.,1000.]
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,2)
ndata= n_elements(indx)
ydata= reform(V[indx],1,ndata)
ycovar= reform(sV[indx]^2D0,1,1,ndata)
projected_gauss_mixtures_c, 2, ydata, ycovar, amp, xmean, xcovar, $
  logfile='tmp'

print, xmean, amp
projection= dblarr(1,1,ndata)+1D0
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,1,2)
assign_clump_members, logpost, ydata, ycovar, projection, xmean, xcovar, amp

IF keyword_set(hardcut) THEN BEGIN
    outer= where(logpost[1,*] GT alog(0.7))
    cfeindx= where(cfe[indx] GT 0.9 and logpost[1,*] GT alog(0.7)) 
    cfeOuter= double(n_elements(cfeindx))/n_elements(outer)
    fehOuter= mean(feh[indx[outer]])
    err= jackknife_frac_hardcut(cfe[indx],cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='red'
ENDIF ELSE BEGIN
    cfeindx= where(cfe[indx] GT 0.9)
    cfeOuter= total(exp(logpost[1,cfeindx]))/total(exp(logpost[1,*]))
    fehOuter= total(exp(logpost[1,*])*feh[indx])/total(exp(logpost[1,*]))
    err= jackknife_frac_weights(cfe[indx],exp(logpost[1,*]),cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='red'
ENDELSE

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='red'
IF keyword_set(hardcut) THEN BEGIN
    xs= feh[indx[outer]]
    predictedCfe= total((coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))/n_elements(outer)
ENDIF ELSE BEGIN
    xs= feh[indx]
    predictedCfe= total(exp(logpost[1,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))
    predictedCfe/= total(exp(logpost[1,*]))
ENDELSE
print, predictedCfe
djs_oplot, [fehOuter], [predictedCfe], psym=4, color='red', symsize=2


IF keyword_set(hardcut) THEN BEGIN
    outer= where(logpost[0,*] GT alog(0.7))
    cfeindx= where(cfe[indx] GT 0.9 and logpost[0,*] GT alog(0.7)) 
    cfeOuter= double(n_elements(cfeindx))/n_elements(outer)
    fehOuter= mean(feh[indx[outer]])
    err= jackknife_frac_hardcut(cfe[indx],cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='green'
ENDIF ELSE BEGIN
    cfeindx= where(cfe[indx] GT 0.9)
    cfeOuter= total(exp(logpost[0,cfeindx]))/total(exp(logpost[0,*]))
    fehOuter= total(exp(logpost[0,*])*feh[indx])/total(exp(logpost[0,*]))
    err= jackknife_frac_weights(cfe[indx],exp(logpost[0,*]),cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='green'
ENDELSE

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='green'
IF keyword_set(hardcut) THEN BEGIN
    xs= feh[indx[outer]]
    predictedCfe= total((coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))/n_elements(outer)
ENDIF ELSE BEGIN
    predictedCfe= total(exp(logpost[0,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))
    predictedCfe/= total(exp(logpost[0,*]))
ENDELSE
print, predictedCfe
djs_oplot, [fehOuter], [predictedCfe], psym=4, color='green', symsize=2


indx= where(feh LT -1.5 and feh GT -2.)
amp= [0.5,0.5]
xmean= [0.,-100.]
xcovar=[1000.,1000.]
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,2)
ndata= n_elements(indx)
ydata= reform(V[indx],1,ndata)
ycovar= reform(sV[indx]^2D0,1,1,ndata)
projected_gauss_mixtures_c, 2, ydata, ycovar, amp, xmean, xcovar, $
  logfile='tmp'

print, xmean, amp
projection= dblarr(1,1,ndata)+1D0
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,1,2)
assign_clump_members, logpost, ydata, ycovar, projection, xmean, xcovar, amp

IF keyword_set(hardcut) THEN BEGIN
    outer= where(logpost[1,*] GT alog(0.7))
    cfeindx= where(cfe[indx] GT 0.9 and logpost[1,*] GT alog(0.7)) 
    cfeOuter= double(n_elements(cfeindx))/n_elements(outer)
    fehOuter= mean(feh[indx[outer]])
    err= jackknife_frac_hardcut(cfe[indx],cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='red'
ENDIF ELSE BEGIN
    cfeindx= where(cfe[indx] GT 0.9)
    cfeOuter= total(exp(logpost[1,cfeindx]))/total(exp(logpost[1,*]))
    fehOuter= total(exp(logpost[1,*])*feh[indx])/total(exp(logpost[1,*]))
    err= jackknife_frac_weights(cfe[indx],exp(logpost[1,*]),cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='red'
ENDELSE


print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='red'
IF keyword_set(hardcut) THEN BEGIN
    xs= feh[indx[outer]]
    predictedCfe= total((coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))/n_elements(outer)
ENDIF ELSE BEGIN
    predictedCfe= total(exp(logpost[1,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))
    predictedCfe/= total(exp(logpost[1,*]))
ENDELSE
print, predictedCfe
djs_oplot, [fehOuter], [predictedCfe], psym=4, color='red', symsize=2

IF keyword_set(hardcut) THEN BEGIN
    outer= where(logpost[0,*] GT alog(0.7))
    cfeindx= where(cfe[indx] GT 0.9 and logpost[0,*] GT alog(0.7)) 
    cfeOuter= double(n_elements(cfeindx))/n_elements(outer)
    fehOuter= mean(feh[indx[outer]])
    err= jackknife_frac_hardcut(cfe[indx],cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='red'
ENDIF ELSE BEGIN
    cfeindx= where(cfe[indx] GT 0.9)
    cfeOuter= total(exp(logpost[0,cfeindx]))/total(exp(logpost[0,*]))
    fehOuter= total(exp(logpost[0,*])*feh[indx])/total(exp(logpost[0,*]))
    err= jackknife_frac_weights(cfe[indx],exp(logpost[0,*]),cempcut=0.9)/100.
    djs_oploterr, [fehOuter], [cfeOuter], yerr=[err], psym=5, color='red'
ENDELSE

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='green'
IF keyword_set(hardcut) THEN BEGIN
    xs= feh[indx[outer]]
    predictedCfe= total((coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))/n_elements(outer)
ENDIF ELSE BEGIN
    predictedCfe= total(exp(logpost[0,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.))
    predictedCfe/= total(exp(logpost[0,*]))
ENDELSE
print, predictedCfe
djs_oplot, [fehOuter], [predictedCfe], psym=4, color='green', symsize=2

END
