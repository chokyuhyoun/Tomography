; Generate Sheep-Logan Phantom


cd, '/data/home/chokh/IDL_default/idl_lib/tomography'
readcol, '/data/home/chokh/IDL_default/idl_lib/tomography/Sheep-Logan Phantom.txt', $
         name, cenx, ceny, major, minor, theta, gray, f='a, f, f, f, f, f'

xs = 128
phantom = fltarr(xs, xs)
xp = (findgen(xs)-0.5*(xs-1))/(0.5*(xs-1))
xxp = rebin(xp, xs, xs)
yyp = rebin(transpose(xp), xs, xs)
inside = fltarr(xs, xs, n_elements(name))
for i=0, n_elements(name)-1 do begin
  rad = theta[i]*!dtor
  dum = ((xxp-cenx[i])*cos(rad)+(yyp-ceny[i])*sin(rad))^2./major[i]^2.$
       +((xxp-cenx[i])*sin(rad)-(yyp-ceny[i])*cos(rad))^2./minor[i]^2.
  dum1 = dum
  dum1[where(dum gt 1.)] = 0.
  dum1[where(dum le 1.)] = 1.
  inside[*, *, i] = dum1
  dum2 = fltarr(xs, xs)
  dum2[where(dum1)] = gray[i]
  phantom = phantom+dum2
;  stop
endfor
save, phantom, filename='SL Phantom.sav'


end