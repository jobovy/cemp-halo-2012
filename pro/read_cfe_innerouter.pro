;+
;
;-
PRO READ_CFE_INNEROUTER, vphi, cfe, p_in, p_out, datafilename=datafilename
IF ~keyword_set(datafilename) THEN datafilename='../data/Test_for_Jo.txt'
OPENR, lun, datafilename, /GET_LUN
hdr= ""
readf, lun, hdr

vphi= dblarr(127)
cfe= dblarr(127)
p_in= dblarr(127)
p_out= dblarr(127)

FOR ii=0L, 126 DO BEGIN
    readf, lun, v,c,pi,po
    vphi[ii]= v
    cfe[ii]= c
    p_in[ii]=pi
    p_out[ii]=po
ENDFOR

END
