;+
;
;-
PRO READ_ERRORDATA, vphi, evphi, datafilename=datafilename, $
                         nhdr=nhdr, ndata=ndata
IF ~keyword_set(datafilename) THEN datafilename='../data/input_dec_ext_d_jo.txt'
IF ~keyword_set(nhdr) THEN nhdr= 2
IF ~keyword_set(ndata) THEN ndata= 344
OPENR, lun, datafilename, /GET_LUN
hdr= ""
FOR ii=0L, nhdr-1 DO readf, lun, hdr

vphi= dblarr(ndata)
evphi= dblarr(ndata)

FOR ii=0L, ndata-1 DO BEGIN
    readf, lun, name, feh, vp, evp, cfe
    vphi[ii]= vp
    evphi[ii]= evp
ENDFOR

END
