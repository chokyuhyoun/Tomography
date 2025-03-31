path = '/data/home/chokh/tomography-1902'
sav_path = path+'/sav00'
file_mkdir, sav_path

cd, path
restore, 'dem_result_5.5_6.5.sav'
restore, 'sdo_indices.sav'
lgt0 = lgtmin+findgen(nlgt)*dlgt
t0=julday(2, 15, 2019, 00, 00, 00)
obs_t = anytim2jul(oindex[*, 0].t_obs)
delt = obs_t-t0
sz=128

;read_sdo, '/data/home/chokh/tomography-1902/sdo_data/094/AIA20190220_140011_0094.fits', $
;          nindex
rot_co=[14.713, -2.396, -1.787]  ;; differential rotation coeff.
xc = oindex[0].crpix1*sz/oindex[0].naxis1
yc = oindex[0].crpix2*sz/oindex[0].naxis2
r_sun=oindex[0].r_sun*sz/oindex[0].naxis1

del_s = oindex[0].naxis1/sz*oindex[0].cdelt1*725d5  ;; 1.4d9 cm

;; tomography
lat=asin((findgen(sz)-yc)/r_sun)
lat[where(~finite(lat) and (findgen(sz) lt xc))]=-0.5*!dpi
lat[where(~finite(lat) and (findgen(sz) gt xc))]=0.5*!dpi
period=360./(rot_co[0]+rot_co[1]*sin(lat)^2.+rot_co[2]*sin(lat)^4.)

cd, sav_path
if 0 then begin
  split_for, 0, nlgt-1, commands=[$   ; temperature
    'ti = systime(/sec)', $
    'box=fltarr(sz, sz, sz)', $
    'print, lgt0[i]', $
    'for j=0l, sz-1 do begin', $  ; latitude
    '  theta=delt/period[j]*360d0', $
    '  box[*, *, j]=solar_tomography(reform(result[*, j, *, i]), $', $
    '               theta, r_sun=r_sun, lat=lat[j])', $
    'endfor', $
    'print, (systime(/sec)-ti)/6d1', $
    'save, box, filename="box_"+string(lgt0[i], f="(f4.2)")+".sav"'], $
    varnames=['lgt0', 'delt', 'period', 'sz', 'r_sun', 'lat', 'result', 'del_s'], $
    ctvariable_name='i', nsplit=12
endif

if 0 then begin
  for i=0, nlgt-1 do begin
    i = 5
    ti = systime(/sec)
    box=fltarr(sz, sz, sz)
    print, lgt0[i]
    for j=0l, sz-1 do begin
      theta=delt/period[j]*360d0
      box[*, *, j]=solar_tomography(reform(result[*, j, *, i]*1d26), $
                   theta, r_sun=r_sun, lat=lat[j])
    endfor
    print, (systime(/sec)-ti)
    stop
  endfor
endif

sav_file = file_search('box_?.??.sav')
nlgt = n_elements(sav_file)
emiss = fltarr(sz, sz, sz, nlgt)
for i=0, nlgt-1 do begin
  restore, sav_file[i]
  box[where(box lt 0.)] = 0.
  emiss[*, *, *, i] = box
endfor
;stop
;nlgt = 11
emiss = emiss[*, *, *, 0:nlgt-1]
 
n_e = sqrt(abs(emiss*1d26/del_s))
n_e[where(emiss lt 0.)] = 0.
tot_e = sqrt(total(emiss*1d26/del_s, 4, /nan))
tbin1 = 10.^(lgtmin+findgen(nlgt)*dlgt)
tbin2 = 10.^(lgtmin+(findgen(nlgt)+1.)*dlgt)
tbin = 0.5*(tbin1+tbin2)
tbinn = rebin(reform(tbin, 1, 1, 1, nlgt), sz, sz, sz, nlgt)
avg_te = total(tbinn*n_e, 4, /nan)/total(n_e, 4, /nan)
avg_te2 = total(tbinn*emiss, 4, /nan)/total(emiss, 4, /nan)

cd, '..'
save, emiss, n_e, tot_e, avg_te, avg_te2, lgtmin, nlgt, dlgt, r_sun, $
      sz, delt, t0, del_s, $
      filename='tdem_result_5.5_6.5.sav'

end