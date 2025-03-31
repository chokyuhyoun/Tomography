path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path

restore, '3d_sinogram_fan_deg7.sav'
sav_file = '3d_reconst_fan_deg7.sav' 
;; head3d, sinogram, theta, ntheta, dtheta, r_sun, tilt_ang
sz = (size(sinogram))[1]
xc = 0.5*(sz-1.)
lat = asin((findgen(sz)-xc)/r_sun)
lat[where(~finite(lat) and (findgen(sz) lt xc))] = -0.5*!dpi
lat[where(~finite(lat) and (findgen(sz) gt xc))] = 0.5*!dpi

box = fltarr(sz, sz, sz)
for i=0, sz-1 do begin
  print, i
  box[*, *, i] = solar_tomography(reform(sinogram[*, i, *]), $
                                  theta, r_sun=r_sun, lat=lat[i])
endfor
save, sz, head3d, box, r_sun, tilt_ang, filename=sav_file

end