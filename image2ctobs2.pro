path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path
data1 =read_png('KakaoTalk_20170110_124327414.png')
data = fltarr(5d2, 5d2)+255
data[*, 0:-2] = data1[0, *, *] 
data[0:200, 400:499] = 255
data[480:*, *] = 255
data = shift(data, [0, 14])
;data = 255-data
sz = 256
r_sun = 44
data=float(congrid(data, sz, sz))
xp=findgen(sz)
yp=findgen(sz)
xc=(sz-1)*0.5
yc=(sz-1)*0.5
xxp=rebin(xp, sz, sz)-xc
yyp=rebin(transpose(yp), sz, sz)-yc
data[where(sqrt(xxp^2.+yyp^2.) le r_sun)]=0.
data[where(sqrt(xxp^2.+yyp^2.) gt sz*0.5)]=0.
dtheta=0.5
ntheta=360./dtheta
theta=dindgen(ntheta)*dtheta-180
sinogram=fltarr(sz, ntheta)
stop
for i=0, ntheta-1 do begin
  rot_data=fiss_embed(data, angle=-theta[i], side=sz)
  dum=rot_data
  dum[where(abs(xxp) lt r_sun and yyp gt 0)]=0.
  sinogram[*, i]=total(dum, 2)
  tvscl, dum
;  stop
endfor
;stop
save, data, sinogram, theta, ntheta, dtheta, r_sun, $
      filename='sinogram2.sav'
end