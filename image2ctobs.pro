path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path

;f=filepath('head.dat', subdir=['examples', 'data'])
;data=read_binary(f, data_dims=[80, 100, 57])
;data=congrid(data[*, *, 30], sz, sz)

;data = read_png('800px-SheppLogan_Phantom.svg.png')
;data = congrid(reform(data[0, *, *]), sz, sz)

;restore, 'SL Phantom.sav'
;sz = (size(phantom))[1]
;data = phantom

restore, 'head3d.sav' ;; head3d, sz
data = reform(head3d[*, *, 64])

xp = findgen(sz)-(sz-1)*0.5
xxp = rebin(xp, sz, sz)
yyp = rebin(transpose(xp), sz, sz)
dist = sqrt(xxp^2.+yyp^2.)
data[where(dist gt sz*0.5)] = 0.
dtheta=0.5
ntheta=180./dtheta
theta=dindgen(ntheta)*dtheta
sinogram=fltarr(sz, ntheta)
for i=0, ntheta-1 do begin
  sinogram[*, i]=total(fiss_embed(data, angle=-theta[i], side=sz), 2)
endfor
save, data, sinogram, theta, ntheta, dtheta, $
      filename='sinogram.sav'
end