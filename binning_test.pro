path = '/hae/homedata/khcho/tomography-1906/sdo_data'
cd, path
f = file_search(path, 'AIA20190529_0000??_*.fits')

img1 = fltarr(128, 128, 6)
img2 = img1
read_sdo, f, index
for i=0, n_elements(f)-1 do begin
  read_sdo, f[i], index1, data1
  img1[*, *, i] = rebin(float(data1), 128, 128)/index1.exptime
  img2[*, *, i] = binning(data1, 32)/index1.exptime/32.^2*2.
endfor

lgTmin = 5.5   ; minimum for lgT axis for inversion
dlgT   = 0.1   ; width of lgT bin
nlgT   = 13    ; number of lgT bins (lgT = [5.5, 6.7]
aia_sparse_em_init, timedepend = index[0].date_obs, /evenorm, $
  use_lgtaxis=findgen(nlgT)*dlgT+lgTmin
lgtaxis = aia_sparse_em_lgtaxis()
exptimestr = '['+strjoin(string(index.exptime,format='(F8.5)'),',')+']'
tolfunc = 'aia_bp_estimate_error(y*'+exptimestr+$
  ', [94,131,171,193,211,335], num_images='+$
  strtrim(string(1^2,format='(I4)'),2)+')/'+exptimestr
aia_sparse_em_solve, img1>0., tolfunc=tolfunc, tolfac=1.4, $
  oem=emcube1, status=status1, coeff=coeff1
aia_sparse_em_solve, img2>0., tolfunc=tolfunc, tolfac=1.4, oem=emcube2, status=status2, coeff=coeff2
end