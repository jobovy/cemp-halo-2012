;+
;
;-
PRO READ_CFE_INNEROUTER, vphi, cfe, p_in, p_out, datafilename=datafilename, $
                         nhdr=nhdr, ndata=ndata
IF ~keyword_set(datafilename) THEN datafilename='../data/Test_for_Jo.txt'
IF ~keyword_set(nhdr) THEN nhdr= 1
IF ~keyword_set(ndata) THEN ndata= 127
OPENR, lun, datafilename, /GET_LUN
hdr= ""
FOR ii=0L, nhdr-1 DO readf, lun, hdr

vphi= dblarr(ndata)
cfe= dblarr(ndata)
p_in= dblarr(ndata)
p_out= dblarr(ndata)

FOR ii=0L, ndata-1 DO BEGIN
    readf, lun, v,c,pi,po
    vphi[ii]= v
    cfe[ii]= c
    p_in[ii]=pi
    p_out[ii]=po
ENDFOR

END
