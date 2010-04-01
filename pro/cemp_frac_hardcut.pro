;+
;   NAME:
;      cemp_frac_hardcut
;   PURPOSE:
;      calculate the fraction of cemp in inner/outer halo as a
;      function of hardcut
;   INPUT:
;      feh - 0 for -1.5, 1 for -2
;      plotfilename
;   OUTPUT:
;      plot
;   HISTORY:
;      2010-04-01 - Written - Bovy (NYU)
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
PRO CEMP_FRAC_HARDCUT, feh=feh, plotfilename=plotfilename
IF ~keyword_set(feh) THEN feh= 0
if ~keyword_set(plotfilename) THEN plotfilename= 'cemp_frac_hardcut_-1.5.ps'
IF feh EQ 0 THEN BEGIN
    READ_CFE_INNEROUTER, vphi, cfe, p_in, p_out, datafilename='../data/Test_for_Jo_-1.5.txt',nhdr=2,ndata=975
    title= '[Fe/H] < -1.5'
    yrange=[0,20]
ENDIF ELSE BEGIN
    READ_CFE_INNEROUTER, vphi, cfe, p_in, p_out, datafilename='../data/Test_for_Jo.txt',nhdr=1,ndata=127
    title= '[Fe/H] < -2'
    yrange=[0,40]
ENDELSE

ncuts= 9
cuts= dindgen(ncuts)/(ncuts-1)*.8
frac_in= dblarr(ncuts)
frac_out= dblarr(ncuts)
efrac_in= dblarr(ncuts)
efrac_out= dblarr(ncuts)
FOR ii=0L, ncuts-1 DO BEGIN
    cfe_in= cfe[where(p_in GT cuts[ii])]
    cfe_out= cfe[where(p_out GT cuts[ii])]
    frac_in[ii]= double(n_elements(cfe_in[where(cfe_in GT 1.)]))/n_elements(cfe_in)*100.
    efrac_in[ii]= jackknife_frac_hardcut(cfe_in)
    frac_out[ii]= double(n_elements(cfe_out[where(cfe_out GT 1.)]))/n_elements(cfe_out)*100.
    efrac_out[ii]= jackknife_frac_hardcut(cfe_out)
ENDFOR

k_print, filename=plotfilename
;;plotting symbol
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill

djs_plot, cuts, frac_in, xtitle='p_{cut}', ytitle='fraction [percent]', xrange=[0,1.], yrange=yrange, color=djs_icolor('green'), title=title
djs_oploterr, cuts, frac_in, yerr=efrac_in, psym=8, color=djs_icolor('green')
djs_oplot, cuts, frac_out, color=djs_icolor('red')
djs_oploterr, cuts, frac_out, yerr=efrac_out, psym=8, color=djs_icolor('red')
k_end_print
END
