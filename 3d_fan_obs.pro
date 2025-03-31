;; 3d sinogram, parallel beam, solar obs

path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path

sav_file = '3d_sinogram_fan_deg7.sav'
tilt_ang = 7.

restore, 'head3d.sav' ;; head3d, sz
r_sun = 50.76

;sz = 32
;r_sun = 10.

xc = 0.5*(sz-1)
obs_dist = -1.5d8/(725.*0.6*32.)   ;; in 32 binning pixel
;obs_dist = -sz
obs_pos = [xc, obs_dist, xc]

xp = findgen(sz)
yyp = rebin(reform(xp, 1, sz), sz, sz, sz)   ;; los direction
xxp = (rebin(xp, sz, sz, sz)-xc)*(obs_dist-(yyp-xc))/obs_dist+xc
zzp = (rebin(reform(xp, 1, 1, sz), sz, sz, sz)-xc)*(obs_dist-(yyp-xc))/obs_dist+xc

xxp2 = rebin(xp, sz, sz)
zzp2 = rebin(transpose(xp), sz, sz)
dist_cen_plane = sqrt((xxp2-xc)^2.+(zzp2-xc)^2.)
cos_plane = abs(obs_dist)/sqrt(dist_cen_plane^2.+obs_dist^2.)
dist_cen_los = dist_cen_plane*cos_plane
blocked = where(dist_cen_los lt r_sun)
 
touch_sun = fltarr(sz, sz)+sz-1 
fan_filter = fltarr(sz, sz, sz)+1.
for i=0, n_elements(blocked)-1 do begin
  los_arr = array_indices(xxp2, blocked[i])
  dist_cen_line = sqrt((xxp[los_arr[0], *, los_arr[1]]-xc)^2.+ $
                       (yyp[los_arr[0], *, los_arr[1]]-xc)^2.+ $
                       (zzp[los_arr[0], *, los_arr[1]]-xc)^2.)
  touch = min(where(dist_cen_line lt r_sun))
  touch_sun[los_arr[0], los_arr[1]] = touch
  fan_filter[los_arr[0], touch:*, los_arr[1]] = 0. 
endfor

dtheta=0.5
ntheta=360./dtheta
theta=(dindgen(ntheta)-0.5*ntheta)*dtheta

sinogram = fltarr(sz, sz, ntheta)
for i=0, ntheta-1 do begin
;    i = 90
  print, i
  rot_ang = [tilt_ang, 0, -theta[i]]
  t3d, /reset, rotate=rot_ang
  pos=!p.t[0:2, 0:2]##[[xxp[*]-xc], [yyp[*]-xc], [zzp[*]-xc]]
  posx=pos[*, 0]+xc
  posy=pos[*, 1]+xc
  posz=pos[*, 2]+xc
  temp_cube=interpolate(head3d, posx, posy, posz, missing=0.)
  t3d, /reset
  temp_cube = reform(temporary(temp_cube), sz, sz, sz)
  temp_cube1 = temp_cube*fan_filter
  sinogram[*, *, i] = total(temp_cube1, 2)
  tvscl, temp_cube[64, *, *]
;  stop
endfor
sinogram = temporary(sinogram)/rebin(cos_plane, sz, sz, ntheta) ;; los facter 
;stop
save, head3d, sinogram, theta, ntheta, dtheta, r_sun, tilt_ang, $
  filename=sav_file
end