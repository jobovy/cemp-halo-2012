PRO MEAN_ERRS, njack=njack, datafilename=datafilename
;IF ~keyword_set(datafilename) THEN datafilename='../data/input_dec_ext_d_jo.txt'
IF ~keyword_set(datafilename) THEN datafilename='../data/input_dec_ext_d.txt'
;;Read data
readcol, datafilename, blah, feh, V, sV, cfe, format='A,D,D,D,D'

;;Run XD
ngauss= 2
amp= [.5,.5]
xmean= [-120., 10.]
xmean= reform(xmean,1,ngauss)
xcovar= [32000.,15000.]
xcovar= reform(xcovar,1,1,ngauss)

ndata= n_elements(V)
ydata= reform(V,1,ndata)
ycovar= reform(sV^2D0,1,1,ndata)
projected_gauss_mixtures_c, 2, ydata, ycovar, amp, xmean, xcovar, avgloglikedata=loglike

print, xmean, amp, loglike

IF ~keyword_set(njack) THEN njack= 90
exd_jackerr, 2, ydata, ycovar, amp, xmean, xcovar, amperr, xmeanerr, xcovarerr, $
  /quiet, njack=njack

print, njack
print, amp, amperr
print, xmean[0], xmeanerr[0]
print, xmean[1], xmeanerr[1]
print, sqrt(xcovar[0]), sqrt(xcovarerr[0])
print, sqrt(xcovar[1]), sqrt(xcovarerr[1])

END
