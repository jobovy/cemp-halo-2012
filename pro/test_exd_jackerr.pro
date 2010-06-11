;+
;   NAME:
;      test_exd_jackerr
;   PURPOSE:
;      test the exd_jackerr.pro code
;   INPUT:
;   OUTPUT:
;   HISTORY:
;      2010-06-11 - Written - Bovy (NYU)
;-
PRO TEST_EXD_JACKERR
datafilename= '../data/File_Test_Carbon.txt'
readcol, datafilename, blah, feh, V, sV, cfe, format='A,D,D,D,D'

;indx= where(feh LT -1.5 and feh GT -2.)
indx= where(feh LT 2. and feh GT -2.5)
amp= [0.5,0.5]
xmean= [0.,-100.]
xcovar=[1000.,1000.]
xmean= reform(xmean,1,2)
xcovar= reform(xcovar,1,2)
ndata= n_elements(indx)
ydata= reform(V[indx],1,ndata)
ycovar= reform(sV[indx]^2D0,1,1,ndata)
projected_gauss_mixtures_c, 2, ydata, ycovar, amp, xmean, xcovar

print, xmean, amp

exd_jackerr, 2, ydata, ycovar, amp, xmean, xcovar, amperr, xmeanerr, xcovarerr, $
  /quiet

print, amp, amperr
print, xmean, xmeanerr
print, xcovar, xcovarerr

END
