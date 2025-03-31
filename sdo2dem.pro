path = '/data/home/chokh/tomography-1902'
cd, path
waves = [94, 131, 171, 193, 211, 335]

restore, '/data/home/chokh/tomography-1902/sdo_data/sdo_indices.sav'
restore, '/data/home/chokh/tomography-1902/sdo_data/sdo_img_all.sav'
sz = (size(img))[1]
nt = (size(img))[3]

;stop


lgTmin = 5.5   ; minimum for lgT axis for inversion
dlgT   = 0.1   ; width of lgT bin
nlgT   = 10    ; number of lgT bins (lgT = [5.5, 7.0]
aia_sparse_em_init, timedepend = oindex[0].date_obs, /evenorm, $
                    use_lgtaxis=findgen(nlgT)*dlgT+lgTmin
;stop
lgtaxis = aia_sparse_em_lgtaxis()
result = fltarr(sz, sz, nt, nlgt)
for i=0, nt-1 do begin
  print, i, nt
  exptimestr = '['+strjoin(string(oindex[i, *].exptime,format='(F8.5)'),',')+']'
  tolfunc = 'aia_bp_estimate_error(y*'+exptimestr+$
            ', [94,131,171,193,211,335], num_images='+$
            strtrim(string(1^2,format='(I4)'),2)+')/'+exptimestr
  aia_sparse_em_solve, reform(img[*, *, i, *])>0., tolfunc=tolfunc, tolfac=1.4, $
                       oem=emcube, status=status, coeff=coeff
;  stop
  result[*, *, i, *] = emcube
endfor
save, result, lgtmin, dlgt, nlgt, sz, nt, $
      filename='dem_result_5.5_6.5.sav'
  
end