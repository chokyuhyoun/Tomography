path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path

sz = 128
r_sun = 20.
;r_sun = 50.76


;f=filepath('head.dat', subdir=['examples', 'data'])
;data=read_binary(f, data_dims=[80, 100, 57])
;data=float(congrid(data[*, *, 30], sz, sz))

;data = read_png('800px-SheppLogan_Phantom.svg.png')
;data = congrid(reform(data[0, *, *]), sz, sz)

;restore, 'SL Phantom.sav'
;data = phantom

restore, 'head3d.sav' ;; head3d, sz
data = reform(head3d[*, *, 64])

xp=findgen(sz)
yp=findgen(sz)
xc=(sz-1)*0.5
yc=(sz-1)*0.5
xxp=rebin(xp, sz, sz)-xc
yyp=rebin(transpose(yp), sz, sz)-yc
data[where(sqrt(xxp^2.+yyp^2.) le r_sun)]=0.
data[where(sqrt(xxp^2.+yyp^2.) gt sz)]=0.
dtheta=0.5
ntheta=360./dtheta
theta=(dindgen(ntheta)-0.5*ntheta)*dtheta
sinogram=fltarr(sz, ntheta)
for i=0, ntheta-1 do begin
  rot_data=fiss_embed(data, angle=-theta[i], side=sz)
  dum=rot_data
  dum[where(abs(xxp) lt r_sun and yyp lt 0)]=0.
  sinogram[*, i]=total(dum, 2)
  tvscl, dum
;  stop
endfor
;stop
save, data, sinogram, theta, ntheta, dtheta, r_sun, $
      filename='sinogram1.sav'
end