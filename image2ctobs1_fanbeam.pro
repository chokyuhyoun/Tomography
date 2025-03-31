path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path

sz = 128
r_sun = 20
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

xc=(sz-1)*0.5
yc=(sz-1)*0.5
xp=findgen(sz)-xc
yp=findgen(sz)-yc
xxp=rebin(xp, sz, sz)
yyp=rebin(transpose(yp), sz, sz)
data[where(sqrt(xxp^2.+yyp^2.) le r_sun)]=0.
data[where(sqrt(xxp^2.+yyp^2.) gt sz)]=0.
dtheta=0.5
ntheta=360./dtheta
theta=(dindgen(ntheta)-0.5*ntheta)*dtheta
sinogram=fltarr(sz, ntheta)

xp = findgen(sz)
yp = xp
xp_pos = rebin(xp, sz, sz) ; x position at y = 63.5 
yyp = rebin(transpose(yp), sz, sz)
slope = 1.5d8/((findgen(sz)-xc)*725.*19.2)
slope1 = rebin(slope, sz, sz)
y_fan = yyp
x_fan = (yyp-xc)/slope1+xp_pos

ll = abs((xc-xp)*slope)/sqrt(slope^2.+1.) ;; distance from solar center to los line
sol_surf = -r_sun*sin(acos(ll/r_sun)+0.5*!dpi-abs(atan(slope)))+63.5
scr_filter = fltarr(sz, sz)+1.
for i=0, n_elements(slope)-1 do begin
  if finite(sol_surf[i]) then scr_filter[i, round(sol_surf[i]):*] = 0.
endfor

;; test
if 0 then begin
p41 = plot(r_sun*cos(findgen(361)*!dtor)+xc, $
           r_sun*sin(findgen(361)*!dtor)+yc, $
           xr=[0, sz], yr=[0, sz], aspect_ratio=1)
for i=0, sz-1 do begin
  los_tip = finite(sol_surf[i]) ? sol_surf[i] : sz-1 
  p42 = plot(x_fan[i, 0:los_tip], y_fan[i, 0:los_tip], over=p41)
endfor
endif
     

for i=0, ntheta-1 do begin
  rot_data=fiss_embed(data, angle=-theta[i], side=sz)
  dum = interpolate(rot_data, x_fan, y_fan)*scr_filter
  sinogram[*, i]=total(dum, 2)/abs(cos(atan(slope)-0.5*!dpi))
  tvscl, dum
;  stop
endfor
;stop
save, data, sinogram, theta, ntheta, dtheta, r_sun, $
      filename='sinogram1_fanbeam.sav'
end