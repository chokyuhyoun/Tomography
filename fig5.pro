cd, '/hae/homedata/khcho/tomography-1906'
restore, 'tdem_result_5.5_6.7.sav' ;; tot_e, avg_te
restore, 'sdo_img_all.sav' ;; img

np = 100
rho2 = [1.05:1.25:0.01]
phi2 = 2*!dpi*randomu(seed1, np)
cost = 2.*randomu(seed2, np)-1.
theta2 = acos(cost)
nrho2 = n_elements(rho2)
nphi2 = n_elements(phi2)
ntheta2 = n_elements(theta2)

rho2p = rebin(rho2, nrho2, nphi2, ntheta2)
phi2p = rebin(reform(phi2, 1, nphi2), nrho2, nphi2, ntheta2)
theta2p = rebin(reform(theta2, 1, 1, ntheta2), nrho2, nphi2, ntheta2)

xxp2 = rho2p*r_sun*sin(theta2p)*cos(phi2p)
yyp2 = rho2p*r_sun*sin(theta2p)*sin(phi2p)
zzp2 = rho2p*r_sun*cos(theta2p)    

xxp2 = reform(temporary(xxp2), nrho2, nphi2, ntheta2)+63.5
yyp2 = reform(temporary(yyp2), nrho2, nphi2, ntheta2)+63.5
zzp2 = reform(temporary(zzp2), nrho2, nphi2, ntheta2)+63.5

emis_surf = fltarr(nrho2, nlgt)
for i=0, nlgt-1 do begin
  dum = interpolate(emiss[*, *, *, i], xxp2, yyp2, zzp2, missing=!values.f_nan)
  emis_surf[*, i] = total(total(dum, 3, /nan), 2, /nan) 
endfor
s_tot_e = sqrt(total(emis_surf*1d26/del_s, 2, /nan))

tbin1 = 10.^(lgtmin+findgen(nlgt)*dlgt)
tbin2 = 10.^(lgtmin+(findgen(nlgt)+1.)*dlgt)
tbin = 0.5*(tbin1+tbin2)
tbinn = rebin(reform(tbin, 1, nlgt), nrho2, nlgt)

s_te = total(tbinn*emis_surf, 2, /nan)/total(emis_surf, 2, /nan)
;stop

; 1. total electron density *********************

sz = 128
rot_co=[14.713, -2.396, -1.787]  ;; differential rotation coeff.
lat=asin((findgen(sz)-63.5)/r_sun)
lat[where(~finite(lat) and (findgen(sz) lt 63.5))]=-0.5*!dpi
lat[where(~finite(lat) and (findgen(sz) gt 63.5))]=0.5*!dpi
ang_vel=(rot_co[0]+rot_co[1]*sin(lat)^2.+rot_co[2]*sin(lat)^4.) ; deg/day

ondisk = where(lat gt -0.5*!dpi and lat lt 0.5*!dpi)
surf_img = fltarr(361, n_elements(ondisk))
surf_lat = lat[ondisk]*!radeg
surf_lon = findgen(361)-0.5*360

for i=0, n_elements(ondisk)-1 do begin
  lat1 = lat[ondisk[i]]
  obs_ang = delt*ang_vel[ondisk[i]]
  lat_pos = interpol(findgen(n_elements(delt)), obs_ang, surf_lon)
  surf_img[*, i] = interpol(img[64, ondisk[i], *, 3], $
                            findgen(n_elements(delt)), lat_pos)  
endfor
surf_img = reverse(surf_img, 1) ; in real, past = right

xp = findgen(128)-64
xxp = rebin(xp, 128, 128, 128)
yyp = rebin(reform(xp, 1, 128), 128, 128, 128)
zzp = rebin(reform(xp, 1, 1, 128), 128, 128, 128)
dist = sqrt(xxp^2.+yyp^2.+zzp^2.)

rho = r_sun*(1.05d0+0.05*findgen(3))
theta = surf_lat*!dtor
phi = surf_lon*!dtor-0.5*!dpi

nrho = n_elements(rho)
ntheta = n_elements(theta)
nphi = n_elements(phi)

rho1 = rebin(rho, nrho, nphi, ntheta)
phi1 = rebin(reform(phi, 1, nphi), nrho, nphi, ntheta)
theta1 = rebin(reform(theta, 1, 1, ntheta), nrho, nphi, ntheta)

xxp1 = rho1*cos(theta1)*cos(phi1) + 63.5
yyp1 = rho1*cos(theta1)*sin(phi1) + 63.5
zzp1 = rho1*sin(theta1) + 63.5

maps = interpolate(tot_e, xxp1, yyp1, zzp1)

w01 = window(dim=[8d2, 8.5d2])

;stop
aia_lct, rr, gg, bb, wave=193
ct193 = transpose([[rr], [gg], [bb]])
im01 = image(surf_img, surf_lon, surf_lat, axis=2, $ 
                /current, pos=[90, 655, 380, 805], /dev, $
                rgb_table=ct193, min=20, max=4d2, xr=[-180, 180], $
                xtickinterval=90, ytickinterval = 30, $
                xtickformat='(i0)', xtickdir=1, ytickdir=1, $
                font_size=12, font_style=1, font_name='malgun gothic', $
                xminor=5, yminor=1, xthick=1.5, ythick=1.5, $
                title='!a(a) SDO/AIA 193 $\AA$ Carrington Map', $
                xtitle='longitude (deg)', ytitle='latitude (deg)')

;stop                
;t01 = text(im01.pos[0]+0.03, im01.pos[1]+0.02, 'AIA 193 !z(00c5)', $
;           /normal, font_size=12, font_style=1, font_name='malgun gothic', $
;           color='white')                
im02 = objarr(nrho)
t02 = objarr(nrho)
for i=0, nrho-1 do begin
  prepos = (i eq 0) ? im01.pos-[0, 0.10, 0, 0.10] : im02[i-1].pos  
  im02[i] = image(alog10(reform(maps[i, *, *])>1.), surf_lon, surf_lat, $
                     /current, axis=2, pos=prepos-[0, 0.18, 0, 0.18], $
                     rgb_table=1, min=7.5, max=8.5, $
                     xtickinterval=90, ytickinterval=30, $
                     xminor=5, yminor=1, xthick=1.5, ythick=1.5, $
                     xtickformat=(i eq nrho-1) ? '(i0)' : '(a1)', $
                     title=(i eq 0) ? '!a(c) Electron Number Density Map' : '', $
                     xtitle=(i eq nrho-1) ? 'longitude (deg)' : '', $
                     ytitle=(i eq nrho-1) ? 'latitude (deg)' : '', $
                     xtickdir=1, ytickdir=1, xr=[-180, 180], $
                     font_size=12, font_style=1, font_name='malgun gothic')
  t02[i] = text(im02[i].pos[0]+0.03, im02[i].pos[1]+0.02, $
                '$\rho$ = '+string(rho[i]/r_sun, f='(f4.2)')+' R$_\odot$', $
                /normal, color='white', $
                font_size=12, font_style=1, font_name='malgun gothic')
endfor
cb1 = colorbar(target=im02[i-1], /relative, pos=[0, -0.5, 1, -0.45], $
               orient=0, textpos=0, thick=1.5, $
               title='log N$_e$ (cm$^{-3}$)', /border, $
               font_style=1, font_size=11, font_name='malgun gothic', $
               ticklen=0.5)

t021 = text(-70, 30, 'A', target=im02[0], /data, color='white', $
            align=0.5, vertical_align=0.5, $
            font_size=10, font_style=0, font_name='malgun gothic') 
t022 = text(75, -30, 'B', target=im02[0], /data, color='white', $
            align=0.5, vertical_align=0.5, $
            font_size=10, font_style=0, font_name='malgun gothic')
t023 = text(110, 30, 'C', target=im02[0], /data, color='white', $
            align=0.5, vertical_align=0.5, $            
            font_size=10, font_style=0, font_name='malgun gothic')


; 2. Electron Temperature *********************
maps = interpolate(avg_te, xxp1, yyp1, zzp1)

for i=0, n_elements(ondisk)-1 do begin
  lat1 = lat[ondisk[i]]
  obs_ang = delt*ang_vel[ondisk[i]]
  lat_pos = interpol(findgen(n_elements(delt)), obs_ang, surf_lon)
  surf_img[*, i] = interpol(img[64, ondisk[i], *, 2], $
    findgen(n_elements(delt)), lat_pos)
endfor
surf_img = reverse(surf_img, 1) ; in real, past = right

 ;stop
aia_lct, rr, gg, bb, wave=171
ct193 = transpose([[rr], [gg], [bb]])
im11 = image(surf_img, surf_lon, surf_lat, axis=2, $
   /current, pos=[490, 655, 780, 805], /dev, $
   rgb_table=ct193, min=0, max=3d2, xr=[-180, 180], $
   xtickinterval=90, ytickinterval = 30, $
   xtickformat='(i0)', xtickdir=1, ytickdir=1, $
   font_size=12, font_style=1, font_name='malgun gothic', $
   xminor=5, yminor=1, xthick=1.5, ythick=1.5, $
   title='!a(b) SDO/AIA 171 $\AA$ Carrington Map', $
   xtitle='longitude (deg)', ytitle='latitude (deg)')
;t11 = text(im11.pos[0]+0.03, im11.pos[1]+0.02, 'AIA 171 !z(00c5)', $
;   /normal, font_size=12, font_style=1, font_name='malgun gothic', $
;   color='white')
im12 = objarr(nrho)
t12 = objarr(nrho)
for i=0, nrho-1 do begin
   prepos = (i eq 0) ? im11.pos-[0, 0.10, 0, 0.10] : im12[i-1].pos
   im12[i] = image(alog10(reform(maps[i, *, *])), surf_lon, surf_lat, $
     /current, axis=2, pos=prepos-[0, 0.18, 0, 0.18], $
     rgb_table=colortable(74, /rev), min=5.9d0, max=6.4d0, $
     xtickinterval=90, ytickinterval = 30, $
     xminor=5, yminor=1, xthick=1.5, ythick=1.5, $
     xtickformat=(i eq nrho-1) ? '(i0)' : '(a1)', $
     title=(i eq 0) ? '!a(d) Electron Temperature Map' : '', $
     xtitle=(i eq nrho-1) ? 'longitude (deg)' : '', $
     ytitle=(i eq nrho-1) ? 'latitude (deg)' : '', $
     xtickdir=1, ytickdir=1, xr=[-180, 180], $
     font_size=12, font_style=1, font_name='malgun gothic')
   t12[i] = text(im12[i].pos[0]+0.03, im12[i].pos[1]+0.02, $
     '$\rho$ = '+string(rho[i]/r_sun, f='(f4.2)')+' R$_\odot$', $
     /normal, color='black', $
     font_size=12, font_style=1, font_name='malgun gothic')
;   t13 = text_bg(t12[i])  
endfor
cb2 = colorbar(target=im12[i-1], /relative, pos=[0, -0.5, 1, -0.45], $
   orient=0, textpos=0, thick=1.5, $
   title='log $\langle T_e \rangle$', /border, $
   font_style=1, font_size=11, font_name='malgun gothic', $
   ticklen=0.5, minor=3)

;stop
w01.save, 'fig5.pdf', page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2, /bitmap
end
