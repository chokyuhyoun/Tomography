cd, '/data/home/chokh/tomography-1902'
restore, 'tdem_result_5.5_6.5.sav' ;; tot_e, avg_te
;restore, 'sdo_img_all.sav' ;; img

; 1. total electron density *********************

sz = 128
lat=asin((findgen(sz)-63.5)/r_sun)
lat[where(~finite(lat) and (findgen(sz) lt 63.5))]=-0.5*!dpi
lat[where(~finite(lat) and (findgen(sz) gt 63.5))]=0.5*!dpi

ondisk = where(lat gt -0.5*!dpi and lat lt 0.5*!dpi)
surf_img = fltarr(361, n_elements(ondisk))
surf_lat = lat[ondisk]*!radeg
surf_lon = findgen(361)-0.5*360

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

tbin = lgtmin+findgen(nlgt)*dlgt
ts = [5.6, 5.8, 6.0, 6.2]
nt = n_elements(ts)
ts_arr = arr_eq(tbin, ts)
te = ts+0.2
te_arr = arr_eq(tbin, te)
maps = fltarr(nrho, nphi, ntheta, nt) 
for i=0, nt-1 do begin
  dum = total(emiss[*, *, *, ts_arr[i]:te_arr[i]], 4)
  dum1 = interpolate(dum, xxp1, yyp1, zzp1)
  maps[*, *, *, i] = reform(dum1, nrho, nphi, ntheta)
endfor
maps = alog10(temporary(maps)+1d-5)+26.-alog10(del_s)

w01 = window(dim=[9d2, 10.3d2], /buffer)
pos0 = [90, 855, 380, 995]
xs = pos0[2]-pos0[0]
ys = pos0[3]-pos0[1]
ymar = 80           
im01 = objarr(nt, nrho)
t01 = objarr(nt, nrho)
t02 = objarr(nt)
for jy = 0, 1 do begin
  for ix = 0, 1 do begin
    kt = jy*2+ix
    for raw = 0, nrho-1 do begin
      pos = pos0-[0, 1, 0, 1]*ys*raw+[1, 0, 1, 0]*410*ix $
                -([0, 1, 0, 1]*ys*3.+[0, 1, 0, 1]*ymar)*jy
      im01[kt, raw] = image(reform(maps[raw, *, *, kt]), $
        surf_lon, surf_lat, /dev, aspect_ratio=0, $
        /current, axis=2, pos=pos, $
        rgb_table=20, min=13, max=17, $
        xtickinterval=90, ytickinterval=30, $
        xminor=5, yminor=1, xthick=1.5, ythick=1.5, $
        xtickformat=(raw eq nrho-1) ? '(i0)' : '(a1)', $
;        title=(i eq 0) ? '!a(c) Electron Number Density Map' : '', $
        xtitle=(raw eq nrho-1) ? 'longitude (deg)' : '', $
        ytitle=(raw eq nrho-1) ? 'latitude (deg)' : '', $
        xr=[-180, 180], $
        font_size=12, font_style=1, font_name='malgun gothic')
      t01[kt, raw] = text(im01[kt, raw].pos[0]+0.03, im01[kt, raw].pos[1]+0.02, $
        '$\rho$ = '+string(rho[raw]/r_sun, f='(f4.2)')+' R$_\odot$', $
        /normal, color='black', $
        font_size=12, font_style=1, font_name='malgun gothic')
      if raw eq 0 then begin
        t02[kt] = text(im01[kt, raw].pos[0]+0.03, im01[kt, raw].pos[3]-0.03, $
          'log $T_e$ = ['+string(ts[kt], f='(f3.1)')+', '+$
                      string(te[kt], f='(f3.1)')+']', $
          /normal, color='black', $
          font_size=12, font_style=1, font_name='malgun gothic')
      endif
;      stop
    endfor
  endfor
endfor
cb1 = colorbar(target=im01[3, 2], /relative, pos=[1, 0, 1.03, 1], $
  orient=1, textpos=1, thick=1.5, minor=1, $
  title='log N$_e^2$ (cm$^{-6}$)', /border, $
  font_style=1, font_size=11, font_name='malgun gothic', $
  ticklen=0.5)

w01.save, 'fig5.pdf', page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2, /bitmap

end
