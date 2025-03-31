cd, '/data/home/chokh/tomography-1902'
restore, 'tdem_result.sav'  ;; box_all
restore, 'avg_te.sav'       ;; avg_te
restore, 'sdo_img_all.sav' ;; img


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
;theta = [-90 : 90 : 1]*!dtor
;phi = [0 : 360 : 1]*!dtor
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

maps = interpolate(avg_te, xxp1, yyp1, zzp1)

;stop
w01 = window(dim=[6d2, 8d2])
aia_lct, rr, gg, bb, wave=193
ct193 = transpose([[rr], [gg], [bb]])
im01 = image(surf_img, surf_lon, surf_lat, axis=2, $ 
                /current, pos=[60, 600, 550, 750], /dev, $
                rgb_table=ct193, min=0, max=5d2, xr=[-180, 180], $
                xtickinterval=90, ytickinterval = 30, $
                xtickformat='(a1)', xtickdir=1, ytickdir=1, $
                font_size=13, font_style=1, font_name='malgun gothic', $
                xminor=5, yminor=1)
t01 = text(im01.pos[0]+0.03, im01.pos[1]+0.02, 'AIA 193 !z(00c5)', $
           /normal, font_size=15, font_style=1, font_name='malgun gothic', $
           color='white')                
im02 = objarr(nrho)
t02 = objarr(nrho)
for i=0, nrho-1 do begin
  prepos = (i eq 0) ? im01.pos : im02[i-1].pos  
  im02[i] = image(alog10(reform(maps[i, *, *])), surf_lon, surf_lat, $
                     /current, axis=2, pos=prepos-[0, 0.21, 0, 0.21], $
                     rgb_table=colortable(74, /rev), min=6.0d0, max=6.4d0, $
                     xtickinterval=90, ytickinterval = 30, $
                     xminor=5, yminor=1, $
                     xtickformat=(i eq nrho-1) ? '(i0)' : '(a1)', $
                     xtitle=(i eq nrho-1) ? 'longitude (deg)' : '', $
                     ytitle=(i eq nrho-1) ? 'latitude (deg)' : '', $
                     xtickdir=1, ytickdir=1, xr=[-180, 180], $
                     font_size=13, font_style=1, font_name='malgun gothic')
  t02[i] = text(im02[i].pos[0]+0.03, im02[i].pos[1]+0.02, $
                '$\rho$ = '+string(rho[i]/r_sun, f='(f4.2)')+' R$_\odot$', $
                /normal, color='black', $
                font_size=15, font_style=1, font_name='malgun gothic')
endfor
im02[2].axes[3].hide = 0
cb1 = colorbar(target=im02[i-1], /relative, pos=[1, 0, 1.04, 1], $
               orient=1, textpos=1, thick=1, $
               title='log(T$_e$) (K)', /border, $
               font_style=1, font_size=11, font_name='malgun gothic', $
               ticklen=0.5, minor=3)

w01.save, 'te_map.pdf', page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2
end