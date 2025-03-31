;; 3d sinogram, parallel beam, solar obs

path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path

restore, 'head3d.sav' ;; head3d, sz
r_sun = 50.76
xc = 0.5*(sz-1)

xp = findgen(sz)
yp = xp
zp = xp
xxp = rebin(xp, sz, sz, sz)
yyp = rebin(reform(yp, 1, sz), sz, sz, sz)   ;; los direction
zzp = rebin(reform(zp, 1, 1, sz), sz, sz, sz)

sol_int = where((xxp-xc)^2.+(yyp-xc)^2.+(zzp-xc)^2. le r_sun^2.)
backside = where(((xxp-xc)^2.+(zzp-xc)^2.) lt r_sun^2. and (yyp ge xc))

dtheta=0.5
ntheta=360./dtheta
theta=(dindgen(ntheta)-0.5*ntheta)*dtheta
stop
sinogram = fltarr(sz, sz, ntheta)
for i=0, ntheta-1 do begin
;  i = 90
  print, i
  temp_cube = rot_cube(head3d, [0, 0, -theta[i]])
  temp_cube[sol_int] = 0.
  temp_cube[backside] = 0.
;  stop
  sinogram[*, *, i] = total(temp_cube, 2)
endfor

save, head3d, sinogram, theta, ntheta, dtheta, r_sun, sol_int, backside, $
      filename='3d_sinogram_parallel.sav'
end