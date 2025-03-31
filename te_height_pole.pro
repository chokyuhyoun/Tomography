cd, '/data/home/chokh/tomography-1902'
restore, 'avg_te.sav' ;; box_all


sz = 128

xp = findgen(128)-64
xxp = rebin(xp, 128, 128, 128)
yyp = rebin(reform(xp, 1, 128), 128, 128, 128)
zzp = rebin(reform(xp, 1, 1, 128), 128, 128, 128)
dist = sqrt(xxp^2.+yyp^2.+zzp^2.)

rho = r_sun*(1.01d0+0.02*findgen(50))
;theta = [-90 : 90 : 30]*!dtor
theta = [-90, 90.]*!dtor
phi = [0 : 360 : 1]*!dtor

nrho = n_elements(rho)
ntheta = n_elements(theta)
nphi = n_elements(phi)

rho1 = rebin(rho, nrho, nphi, ntheta)
phi1 = rebin(reform(phi, 1, nphi), nrho, nphi, ntheta)
theta1 = rebin(reform(theta, 1, 1, ntheta), nrho, nphi, ntheta)

xxp1 = rho1*cos(theta1)*cos(phi1) + 63.5
yyp1 = rho1*cos(theta1)*sin(phi1) + 63.5
zzp1 = rho1*sin(theta1) + 63.5

maps = interpolate(avg_te, xxp1, yyp1, zzp1, missing=!values.f_nan)

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
avg_te1 = interpolate(avg_te, xxp2, yyp2, zzp2)

avg_te2 = mean(avg_te1, dim=1, /nan)                    

w01 = window(dim=[1d3, 5d2])

p12 = objarr(ntheta)
t12 = objarr(ntheta)
te_pro = mean(maps, dim=2, /nan)
te_pro = mean(te_pro, dim=2, /nan)
;te_pro_std = stddev(maps, dim=2, /nan)

p11 = plot(rho/r_sun, (avg_te2), '-3', /sym_filled, /current, /dev, $
           pos=[150, 70, 550, 470], xthick=2, ythick=2, $
           xminor=3, yminor=3, xticklen=0.04, $
           xr=[1.03, 1.4], yr=10^[5.5, 6.5], $
           xtitle='Distance (R$_\odot$)', $
           ytitle='T$_e$ (K)', $
           font_size=13, font_name='malgun gothic')
t11 = text(1.1, 0.9, /relative, target=p11, $
           'Averaged', font_size=13, font_name='malgun gothic')

ct2 = colortable(39, /transpose)
i=0
  p12[i] = plot(rho/r_sun, (te_pro[*, i]), $
                over=p11, '-', $
                color=ct2[*, 254./ntheta*(i+1)], thick=2, transp=0)
  t12[i] = text(1.1, 0.84-0.06*i, /relative, target=p12[i], $
                'Lat = '+string(theta[ntheta-1-i]*!radeg, f='(i3)')+'$\deg$', $
                color=ct2[*, 254./ntheta*(i+1)], $
                font_size=13, font_name='malgun gothic')
i++


readcol, 'previous/CH_David_1998.txt', r, lgte
p20 = plot(r, 10.^lgte, ' o', over=p12[0])
t20 = legend(position=[1.1, 0.84-0.06*i], /relative, target=p20, $
             color='white', label='David et al. (1998)', $
             font_size=13, font_name='malgun gothic', $
             vertical_align=0.5, horizontal_align=0, $
             sample_width=0) 
i++

readcol, 'previous/CH_Reginald_2009.txt', r, lgte
p21 = plot(r, 10.^lgte, ' X', over=p12[0])
t21 = legend(position=[1.1, 0.84-0.06*i], /relative, target=p21, $
              color='white', label='Reginald et al. (2009)', $
              font_size=13, font_name='malgun gothic', $
              vertical_align=0.5, horizontal_align=0, $
              sample_width=0)
i++

readcol, 'previous/CH_Fisher_1995.txt', r, lgte
p22 = plot(r, 10.^lgte, ' D', over=p12[0])
t22 = legend(position=[1.1, 0.84-0.06*i], /relative, target=p22, $
              color='white', label='Fisher (1995)', $
              font_size=13, font_name='malgun gothic', $
              vertical_align=0.5, horizontal_align=0, $
              sample_width=0)
i++

readcol, 'previous/CH_Ko_1997.txt', r, lgte
p23 = plot(r, 10.^lgte, ' s', over=p12[0])
t23 = legend(position=[1.1, 0.84-0.06*i], /relative, target=p23, $
            color='white', label='Ko et al. (1997)', $
            font_size=13, font_name='malgun gothic', $
            vertical_align=0.5, horizontal_align=0, $
            sample_width=0)
i++

;readcol, 'previous/CH_Marsch_2000.txt', r, lgte
;p23 = plot(r, 10.^lgte, ' tu', over=p12[0])
;t23 = legend(position=[1.1, 0.84-0.06*i], /relative, target=p23, $
;              color='white', label='Marsch et al. (2000)', $
;              font_size=13, font_name='malgun gothic', $
;              vertical_align=0.5, horizontal_align=0, $
;              sample_width=0)
;i++



;stop
w01.save, 'te_height_pole.png', resol=200, /bitmap
end