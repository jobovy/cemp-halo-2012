;+
;   NAME:
;      feh_cemp
;   PURPOSE:
;      derive and fit the cemp vs. feh relation
;   INPUT:
;   OUTPUT:
;   HISTORY:
;      2010-06-01 - Written - Bovy (NYU)
;-
PRO FEH_CEMP
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
xcovar=[100.,100.]
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,2)
ndata= n_elements(indx)
ydata= reform(V[indx],1,ndata)
ycovar= reform(sV[indx]^2D0,1,1,ndata)
projected_gauss_mixtures_c, 2, ydata, ycovar, amp, xmean, xcovar, $
  logfile='tmp'

print, xmean
projection= dblarr(1,1,ndata)+1D0
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,1,2)
assign_clump_members, logpost, ydata, ycovar, projection, xmean, xcovar, amp

cfeindx= where(cfe[indx] GT 0.9)
cfeOuter= total(exp(logpost[1,cfeindx]))/total(exp(logpost[1,*]))
fehOuter= total(exp(logpost[1,*])*feh[indx])/total(exp(logpost[1,*]))

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='red'
xs= feh[indx]
predictedCfe= exp(logpost[1,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.)
predictedCfe/= exp(logpost[1,*])
djs_oplot, [fehOuter], [cfeOuter], psym=1, color='red', symsize=2


cfeindx= where(cfe[indx] GT 0.9)
cfeOuter= total(exp(logpost[0,cfeindx]))/total(exp(logpost[0,*]))
fehOuter= total(exp(logpost[0,*])*feh[indx])/total(exp(logpost[0,*]))

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='green'
predictedCfe= exp(logpost[0,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.)
predictedCfe/= exp(logpost[0,*])
djs_oplot, [fehOuter], [cfeOuter], psym=1, color='green', symsize=2


indx= where(feh LT -1.5 and feh GT -2.)
amp= [0.5,0.5]
xmean= [0.,-100.]
xcovar=[100.,100.]
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,2)
ndata= n_elements(indx)
ydata= reform(V[indx],1,ndata)
ycovar= reform(sV[indx]^2D0,1,1,ndata)
projected_gauss_mixtures_c, 2, ydata, ycovar, amp, xmean, xcovar, $
  logfile='tmp'

print, xmean
projection= dblarr(1,1,ndata)+1D0
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,1,2)
assign_clump_members, logpost, ydata, ycovar, projection, xmean, xcovar, amp

cfeindx= where(cfe[indx] GT 0.9)
cfeOuter= total(exp(logpost[1,cfeindx]))/total(exp(logpost[1,*]))
fehOuter= total(exp(logpost[1,*])*feh[indx])/total(exp(logpost[1,*]))

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='red'
predictedCfe= exp(logpost[1,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.)
predictedCfe/= exp(logpost[1,*])
djs_oplot, [fehOuter], [cfeOuter], psym=1, color='red', symsize=2

cfeindx= where(cfe[indx] GT 0.9)
cfeOuter= total(exp(logpost[0,cfeindx]))/total(exp(logpost[0,*]))
fehOuter= total(exp(logpost[0,*])*feh[indx])/total(exp(logpost[0,*]))

print, fehOuter, cfeOuter
djs_oplot, [fehOuter], [cfeOuter], psym=5, color='green'
predictedCfe= exp(logpost[0,*])*(coeffs[0]+coeffs[1]*xs+coeffs[2]*xs^2.)
predictedCfe/= exp(logpost[0,*])
djs_oplot, [fehOuter], [cfeOuter], psym=1, color='green', symsize=2

END
