cd, '/data/home/chokh/tomography-1902'
restore, 'tdem_result.sav' ;; box_all

; 1. total electron density

sz = 128

xp = findgen(128)-64
xxp = rebin(xp, 128, 128, 128)
yyp = rebin(reform(xp, 1, 128), 128, 128, 128)
zzp = rebin(reform(xp, 1, 1, 128), 128, 128, 128)
dist = sqrt(xxp^2.+yyp^2.+zzp^2.)

rho = r_sun*(1.01d0+0.02*findgen(50))
theta = [-90 : 90 : 30]*!dtor
phi = [0 : 360 : 10]*!dtor

nrho = n_elements(rho)
ntheta = n_elements(theta)
nphi = n_elements(phi)

rho1 = rebin(rho, nrho, nphi, ntheta)
phi1 = rebin(reform(phi, 1, nphi), nrho, nphi, ntheta)
theta1 = rebin(reform(theta, 1, 1, ntheta), nrho, nphi, ntheta)

xxp1 = rho1*cos(theta1)*cos(phi1) + 63.5
yyp1 = rho1*cos(theta1)*sin(phi1) + 63.5
zzp1 = rho1*sin(theta1) + 63.5

tot_e = total(box_all, 4)
maps = interpolate(tot_e, xxp1, yyp1, zzp1, missing=!values.f_nan)

np = 1d4
phi2 = 2*!dpi*randomu(seed1, np)
cost = 2d0*randomu(seed2, np)-1d0
theta2 = acos(cost)

xxp2 = fltarr(np, nrho)
yyp2 = fltarr(np, nrho)
zzp2 = fltarr(np, nrho)

for i=0, nrho-1 do begin
  xxp2[*, i] = rho[i]*sin(theta2)*cos(phi2)+63.5
  yyp2[*, i] = rho[i]*sin(theta2)*sin(phi2)+63.5
  zzp2[*, i] = rho[i]*cos(theta2)+63.5
endfor
avg_ne1 = interpolate(tot_e, xxp2, yyp2, zzp2)

avg_ne = mean(avg_ne1, dim=1, /nan)                    

w01 = window(dim=[1d3, 5d2])
rr = [1:1.8:1d-3]
m_name = ['Newkirk', 'Baumbach', 'Cram', 'Leblanc']
den = fltarr(n_elements(rr), n_elements(m_name))
den[*, 0] = 4.2d4*10.^(4.32/rr) ; Newkirk
den[*, 1] = 1d8*(0.036*rr^(-1.5)+1.55*rr^(-6.)+2.99*rr^(-16.)) ; Baumbach
den[*, 2] = 1.67*10.^(4+4.04/rr) ; Cram
den[*, 3] = 3.3d5*rr^(-2.)+4.1d6*rr^(-4.)+8d7*rr^(-6.) ; Leblanc

m_color = ['red', 'orange', 'green', 'blue']
m_thick = 2
m_linestyle = 0

p01 = objarr(n_elements(m_name))
t01 = objarr(n_elements(m_name))
i=0
p01[0] = plot(rr, den[*, i]*1d-8, /current, pos=[70, 70, 470, 470], /dev, $ 
              /nodata, xthick=2, ythick=2, xticklen=0.04, $
              color=m_color[i], thick=m_thick, linestyle=m_linestyle, $
              xr=[1.03, 1.6], yr=[0, 7], $
              xminor=3, yminor=3, $
              xtitle='Distance (R$_\odot$)', $
              ytitle='Electron Density (10$^8$ cm$^{-3}$)', $
              font_size=15, font_name='malgun gothic')
for i=0, n_elements(m_name)-1 do begin
  p01[i] = plot(rr, den[*, i]*1d-8, over=p01[0], $
                color=m_color[i], thick=m_thick, linestyle=m_linestyle) 
  t01[i] = text(0.65, 0.9-0.06*i, /relative, target=p01[i], $
                m_name[i], color=m_color[i], $
                font_size=13, font_name='malgun gothic')
endfor
p02 = plot(rho/r_sun, avg_ne*1d-8, '-k4', over=p01[0])
t02 = text(0.65, 0.9-0.06*i, /relative, target=p02, $
           'This Work', color='black', $
           font_size=13, font_name='malgun gothic')

p12 = objarr(ntheta)
t12 = objarr(ntheta)
den_pro = mean(maps, dim=2, /nan)

p11 = plot(rho/r_sun, avg_ne*1d-8, '-k4', /current, /dev, $
           pos=[570, 70, 970, 470], xthick=2, ythick=2, $
           xminor=3, yminor=3, xticklen=0.04, $
           xr=[1.03, 1.6], yr=[0, 7], $
           xtitle='Distance (R$_\odot$)', $
           ytitle='Electron Density (10$^8$ cm$^{-3}$)', $
           font_size=13, font_name='malgun gothic')
t11 = text(0.65, 0.9, /relative, target=p11, $
           'Averaged', font_size=13, font_name='malgun gothic')

ct2 = colortable(39, /transpose)
for i=0, ntheta-1 do begin
  p12[i] = plot(rho/r_sun, den_pro[*, i]*1d-8, over=p11, $
                color=ct2[*, 254./ntheta*(i+1)], thick=2, transp=0)
  t12[i] = text(0.65, 0.84-0.06*i, /relative, target=p12[i], $
                'Lat = '+string(theta[ntheta-1-i]*!radeg, f='(i3)')+'$\deg$', $
                color=ct2[*, 254./ntheta*(i+1)], $
                font_size=13, font_name='malgun gothic')
endfor

w01.save, 'e_den_height.png', resol=200, /bitmap
end