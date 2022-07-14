cd, '/hae/homedata/khcho/tomography-1906'
restore, 'tdem_result_5.5_6.7.sav' ;; tot_e, avg_te


w01 = window(dim=[8d2, 8.5d2], background_color='black')

bg01 = image(fltarr(10, 10)+255, pos=[400, 450, 800, 850], /dev, /current, min=0)
bg02 = image(fltarr(10, 10)+255, pos=[0, 0, 800, 450], $
             aspect_ratio=0, /dev, /current, min=0)

;stop
ratio = 0.7*0.5
;pos1 = 0.5*[1.-ratio, 1.-ratio, 1.+ratio, 1.+ratio]-[0.05, 0, 0.05, 0]
pos1 = 0.5*[0.5-ratio, 0.5-ratio, 0.5+ratio, 0.5+ratio]$
      -[0.07, -0.92, 0.07, -0.92]*0.5
disk_in = plot3d(indgen(2), indgen(2), indgen(2), $
                /current, axis=0, clip=0, background_transp=100, $
                xr=64.*[-1, 1], yr=64.*[-1, 1], zr=64.*[-1, 1], $
                aspect_z=1, aspect_ratio=1)
disk_in.pos=pos1

restore, 'sdo_img_all.sav' ;; img

sz = 128
rot_co=[14.713, -2.396, -1.787]  ;; differential rotation coeff.
lat=asin((findgen(sz)-63.5)/r_sun)
lat[where(~finite(lat) and (findgen(sz) lt 63.5))]=-0.5*!dpi
lat[where(~finite(lat) and (findgen(sz) gt 63.5))]=0.5*!dpi
ang_vel=(rot_co[0]+rot_co[1]*sin(lat)^2.+rot_co[2]*sin(lat)^4.) ; deg/day

ondisk = where(lat gt -0.5*!dpi and lat lt 0.5*!dpi)
im_sz = [361, n_elements(ondisk)]
surf_img = fltarr(im_sz[0], im_sz[1])
surf_lat = lat[ondisk]*!radeg
surf_lon = findgen(361)-0.5*360

for i=0, im_sz[1]-1 do begin
  lat1 = lat[ondisk[i]]
  obs_ang =delt*ang_vel[ondisk[i]]
  lat_pos = interpol(findgen(n_elements(delt)), obs_ang, surf_lon)
  surf_img[*, i] = interpol(img[64, ondisk[i], *, 3], $
                            findgen(n_elements(delt)), lat_pos)  
endfor
surf_img = reverse(surf_img, 1) ; in real, past = right

w02 = window(dim=[720, 360], /buffer)
aia_lct, rr, gg, bb, wave=193
im91 = image(surf_img, findgen(im_sz[0])/im_sz[0]*360, $
             findgen(im_sz[1])/im_sz[1]*180.-90., $
             /current, pos=[0, 0, 1, 1], $
             rgb_table=[[rr], [gg], [bb]], min=0, max=4d2)
;stop
for i=0, 330, 30 do begin
  p91 = plot(i*[1, 1], [-90, 90], /data, over=im91, $
             color='light gray', thick=1)
endfor
for i=-60, 60, 30 do begin
  p92 = plot([0, 360], i*[1, 1], /data, over=im91, $
             color='light gray', thick=1)
endfor
surf_data = w02.copywindow(border=0, resolution=300)

oImage = obj_new('IDLgrImage', shift(surf_data, 0, (size(surf_data))[2]*0.25, 0))
mesh_obj, 4, vert, poly, replicate(r_sun, 101, 101)
   
vec = findgen(101)/100.
tex_coord = fltarr(2, 101, 101)
tex_coord[0, *, *] = rebin(vec, 101, 101)
tex_coord[1, *, *] = rebin(transpose(vec), 101, 101)
oSolar = obj_new('IDLgrPolygon', data=vert, polygons=poly, $
                  color=[255, 255, 255], $
                  texture_coord=tex_coord, texture_map=oImage, $
                  /texture_interp)
surf = plot3d([0], [0], [0], over=disk_in, $
              aspect_z=1, sym_object=oSolar, /data, sym_size=ratio)
;stop

arrow_arr = fltarr(3, 5)
arrow_arr[0, *] = [0, 0.5, 0.5, 3, 0]
arrow_arr[2, *] = [-80, -80, 75, 75, 80]
mesh_obj, 6, vts, pol, arrow_arr, p1=100, /closed

vts = transpose(vts)
Oarrow1 = obj_new('IDLgrPolygon', transpose(vts), polygons=pol, $
                  color=[200, 200, 200])
arrow1 = plot3d([0], [0], [0], over=disk_in, clip=0, $
                sym_object=Oarrow1, /data, sym_size=ratio)
t221 = text(0, 0, 85, 'Rotation Axis', /data, target=disk_in, $
            color='white', align=0.5, vertical_align=0.5, clip=0, $
            font_size=12, font_style=1, font_name='malgun gothic', $
            baseline=[1, -1, 0], updir=[0, 0, 1])
    
t3d, /reset, rotate=[90, 0, 0]
vts1 = !p.t[0:2, 0:2]##[[vts[*, 0]], [vts[*, 1]], [vts[*, 2]]]
Oarrow2 = obj_new('IDLgrPolygon', transpose(vts1), polygons=pol, $
                color=[200, 200, 200])
arrow2 = plot3d([0], [0], [0], over=disk_in, clip=0, $
                sym_object=Oarrow2, /data, sym_size=ratio)
t222 = text(0, -80, 0, 'Reference !cTime !cDirection', /data, target=disk_in, $
           color='white', vertical_align=0.5, clip=0, $
           font_size=12, font_style=1, font_name='malgun gothic', $
           baseline=[0, -1, 0], updir=[0, 0, 1])

wid = 10. ; in phi degree
np = 1d2
rho = findgen(25)*0.01+1.
phi = 2*!dpi*randomu(seed1, np)
cost = (1.-cos(wid*!dtor))*randomu(seed2, np)+cos(wid*!dtor)
theta = acos(cost)
nrho = n_elements(rho)

rhop = rebin(rho, nrho, np)
phip = rebin(transpose(phi), nrho, np)
thetap = rebin(transpose(theta), nrho, np)

xxp2 = rhop*r_sun*sin(thetap)*cos(phip)
yyp2 = rhop*r_sun*sin(thetap)*sin(phip)
zzp2 = rhop*r_sun*cos(thetap)    

xxp2 = reform(temporary(xxp2), nrho, np)
yyp2 = reform(temporary(yyp2), nrho, np)
zzp2 = reform(temporary(zzp2), nrho, np)

cone_arr = fltarr(3, 2)
cone_arr[0, *] = [0, r_sun*max(rho)*sin(wid*!dtor)]
cone_arr[2, *] = [0, r_sun*max(rho)]
mesh_obj, 6, vts, pol, cone_arr, p1=100, /closed
vts = transpose(vts)

;; ---------- region D
xrot1 = 110.
yrot1 = 0.
zrot1 = 0.
t3d, /reset, rotate=[xrot1, yrot1, zrot1]  ;; positive = counterclockwise 
vts1 = !p.t[0:2, 0:2]##[[vts[*, 0]], [vts[*, 1]], [vts[*, 2]]]
qs_pos = !p.t[0:2, 0:2]##[[xxp2[*]], [yyp2[*]], [zzp2[*]]]+63.5
;p31 = plot3d(qs_pos[*, 0], qs_pos[*, 1], qs_pos[*, 2], over=disk_in, '.w')
cone_obj1 = obj_new('IDLgrPolygon', transpose(vts1), polygons=pol, $
                    color=[100, 100, 100], alpha_channel=1)
cone1 = plot3d([0], [0], [0], over=disk_in, clip=0, $
               sym_object=cone_obj1, /data, sym_size=ratio)
txtpos = !p.t##[0, 0, 70, 1]
t91 = text(txtpos[0], txtpos[1], txtpos[2], 'D', /data, target=disk_in, $
  color='white', align=0.5, vertical_align=0.5, clip=0, $
  font_size=13, font_style=1, font_name='malgun gothic', $
  baseline=[1, -1, 0], updir=[0, 0, 1])
;stop

;; ----------- region E
xrot2 = 180.
yrot2 = 0. 
zrot2 = 0.
t3d, /reset, rotate=[xrot2, yrot2, zrot2]
ch_pos = !p.t[0:2, 0:2]##[[xxp2[*]], [yyp2[*]], [zzp2[*]]]+63.5
;p32 = plot3d(ch_pos[*, 0], ch_pos[*, 1], ch_pos[*, 2], over=disk_in, '.w')
vts2 = !p.t[0:2, 0:2]##[[vts[*, 0]], [vts[*, 1]], [vts[*, 2]]]
cone_obj2 = obj_new('IDLgrPolygon', transpose(vts2), polygons=pol, $
                    color=[100, 100, 100], alpha_channel=1)
cone2 = plot3d([0], [0], [0], over=disk_in, clip=0, $
                    sym_object=cone_obj2, /data, sym_size=ratio)
txtpos = !p.t##[0, 10, 75, 1]
t92 = text(txtpos[0], txtpos[1], txtpos[2], 'E', /data, target=disk_in, $
            color='white', align=0.5, vertical_align=0.5, clip=0, $
            font_size=13, font_style=1, font_name='malgun gothic', $
            baseline=[1, -1, 0], updir=[0, 0, 1])

disk_in.rotate, -130, /x, /reset
disk_in.rotate, 50, /z



;*********** Density Comparison

rf_size =9
rs_size = 0.5

qs_den1 = interpolate(tot_e, qs_pos[*, 0], qs_pos[*, 1], qs_pos[*, 2])
qs_den = mean(reform(qs_den1, nrho, np), dim=2)
ch_den1 = interpolate(tot_e, ch_pos[*, 0], ch_pos[*, 1], ch_pos[*, 2])
ch_den = mean(reform(ch_den1, nrho, np), dim=2)

rr = [1:1.8:1d-3]
m_name = ['Newkirk', 'Baumbach', 'Cram', 'Guhathakurta']
den = fltarr(n_elements(rr), n_elements(m_name))
den[*, 0] = 4.2d4*10.^(4.32/rr) ; Newkirk
den[*, 1] = 1d8*(0.036*rr^(-1.5)+1.55*rr^(-6.)+2.99*rr^(-16.)) ; Baumbach
den[*, 2] = 1.67*10.^(4+4.04/rr) ; Cram
;den[*, 3] = 3.3d5*rr^(-2.)+4.1d6*rr^(-4.)+8d7*rr^(-6.) ; Leblanc
den[*, 3] = (1736.9*rr^(-13.72)+19.95*rr^(-4.09)+1.316*rr^(-2.))*1d5 ; Guhathakurta 

m_color = ['red', 'orange', 'green', 'blue']
m_thick = 1
m_linestyle = 5

p01 = objarr(n_elements(m_name))
t01 = objarr(n_elements(m_name))
i=0
p19 = plot(indgen(2), current=w01, pos=[480, 505, 770, 795], /dev, $
  /nodata, xthick=1.5, ythick=1.5, xticklen=0.04, $
  color=m_color[i], thick=m_thick, linestyle=m_linestyle, $
  xr=[1.03, 1.25], yr=[0, 2], ylog=0, $
  xminor=4, title='(b) Electron Density', $
  xtitle='Height (R$_\odot$)', $
  ytitle='Electron Density (10$^8$ cm$^{-3}$)', $
  font_size=12, font_style=1, font_name='malgun gothic')
for i=0, n_elements(m_name)-1 do begin
  p01[i] = plot(rr, den[*, i]*1d-8, over=p19, $
    color=m_color[i], thick=m_thick, linestyle=m_linestyle)
  t01[i] = text(0.65, 0.9-0.05*i, /relative, target=p01[i], $
    m_name[i], color=m_color[i], $
    font_size=rf_size, font_name='malgun gothic')
endfor
i++
p02 = plot(rho, qs_den*1d-8, '-k3', over=p01[0])
t02 = text(0.65, 0.9-0.05*i, /relative, target=p02, $
          'Region D', color='black', $
          font_size=rf_size, font_name='malgun gothic')
i++          
p03 = plot(rho, ch_den*1d-8, '-3', color='gray', over=p01[0])
t03 = text(0.65, 0.9-0.05*i, /relative, target=p02, $
          'Region E', color='gray', $
          font_size=rf_size, font_name='malgun gothic')

          
;*********** Temperatuer Comparison D
ct_a = colortable(m_color, ncolors=4, /transp)
i = 0
;nlgt = 11  ;===============================
tbin1 = 10.^(lgtmin+findgen(nlgt)*dlgt)
tbin2 = 10.^(lgtmin+(findgen(nlgt)+1.)*dlgt)
tbin = 0.5*(tbin1+tbin2)
tbinn = rebin(tbin, nlgt, nrho)

qs_emiss = fltarr(nlgt, nrho)
ch_emiss = fltarr(nlgt, nrho)
for ii=0, nlgt-1 do begin
  dum11 = interpolate(emiss[*, *, *, ii], $
                     qs_pos[*, 0], qs_pos[*, 1], qs_pos[*, 2])
  dum12 = reform(dum11, nrho, np) 
  qs_emiss[ii, *] = total(dum12, 2, /nan) 
  
  dum21 = interpolate(emiss[*, *, *, ii], $
                      ch_pos[*, 0], ch_pos[*, 1], ch_pos[*, 2])
  dum22 = reform(dum21, nrho, np) 
  ch_emiss[ii, *] = total(dum22, 2, /nan)
endfor

;qs_emiss = qs_emiss[0:10, *]
;ch_emiss = ch_emiss[0:10, *]
;tbinn = tbinn[0:10, *]
qs_te = total(tbinn*qs_emiss, 1)/total(qs_emiss, 1) 
ch_te = total(tbinn*ch_emiss, 1)/total(ch_emiss, 1)

qs_em = transpose(qs_emiss)*1d26/del_s/np
ch_em = transpose(ch_emiss)*1d26/del_s/np

p29 = image(alog10(qs_em), rho, alog10(tbin1), axis=2, current=w01, $
  pos=[80, 110, 370, 400], /dev, aspect_ratio=0, min=14., max=16, $ 
  xthick=1.5, ythick=1.5, xticklen=0.04, xstyle=1, ystyle=1, $
  xr=[1.03, 1.25], yr=[5.5, 6.7], rgb_table=colortable(0, /rev), $
  xminor=4, yminor=4, title='(c) Electron Temperature in QS',  $
  xtitle='Height (R$_\odot$)', ytitle='log $T_e$', $
  font_size=12, font_style=1, font_name='malgun gothic', transp=20)
cb29 = colorbar(target=p29, position=[0.0, -0.21, 1, -0.18], /relative, $
                border=1, major=3, title='log $n_{e,j}^2$', $
                font_style=1, font_size=11, font_name='malgun gothic', $
                ticklen=0.5, thick=1.5, transp=p29.transp, minor=3)

p24 = plot(rho, alog10(qs_te), '-k3', over=p29)

tpos1 = [0.45, 0.27]
ratio1 = 0.05   
readcol, 'previous/QS_Wheatland_1997.txt', r, lgte1, lgte2
p201 = plot(r, lgte1, ' tu', over=p29, sym_size=rs_size, color=ct_a[*, i])
p202 = plot(r, lgte2, ' o', over=p29, sym_size=rs_size, color=ct_a[*, i])
t20 = text(tpos1[0], tpos1[1]-ratio1*i, 'Wheatland et al. (1997)', /relative, target=p29, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

readcol, 'previous/QS_David_1998.txt', r, lgte
p21 = plot(r, lgte, ' o', over=p29, sym_size=rs_size, color=ct_a[*, i])
t21 = text(tpos1[0], tpos1[1]-ratio1*i, 'David et al. (1998)', /relative, target=p29, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++


readcol, 'previous/QS_Warren_1999.txt', r, lgte1, lgte2
p221 = plot(r, lgte1, ' tu', over=p29, sym_size=rs_size, color=ct_a[*, i])
p222 = plot(r, lgte2, ' o', over=p29, sym_size=rs_size, color=ct_a[*, i])
t21 = text(tpos1[0], tpos1[1]-ratio1*i, 'Warren (1999)', /relative, target=p29, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

readcol, 'previous/QS_Reginald_2009.txt', r, lgte
p23 = plot(r, lgte, ' o', over=p29, sym_size=rs_size, color=ct_a[*, i])
t23 = text(tpos1[0], tpos1[1]-ratio1*i, 'Reginald et al. (2009)', /relative, target=p29, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

t24 = text(tpos1[0], tpos1[1]-ratio1*i, 'Region D', /relative, target=p29, $
           color='black', font_style=1, font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)

  
;*********** Temperatuer Comparison E
ct_a = colortable(['red', 'orange', 'green', 'aqua', 'blue', 'purple'], ncolors=6, /transp)
i = 0
p39 = image(alog10(ch_em), rho, alog10(tbin1), axis=2, current=w01, $
  pos=[480, 110, 770, 400], /dev, aspect_ratio=0, min=13, max=15, $ 
  xthick=1.5, ythick=1.5, xticklen=0.04, xstyle=1, ystyle=1, $
  xr=p29.xr, yr=[5.5, 6.7], rgb_table=colortable(0, /rev), $
  xminor=4, yminor=4, title='(d) Electron Temperature in CH',  $
  xtitle='Height (R$_\odot$)', ytitle='log $T_e$', $
  font_size=12, font_style=1, font_name='malgun gothic', transp=p29.transp)
cb39 = colorbar(target=p39, position=[0.0, -0.22, 1, -0.19], /relative, $
                border=1, major=cb29.major, title='log $n_{e,j}^2$', $
                font_style=1, font_size=11, font_name='malgun gothic', $
                ticklen=0.5, thick=1.5, transp=p39.transp, minor=3)

p36 = plot(rho, alog10(ch_te), '-k3', over=p39)
tpos2 = [0.05, 0.23] 
ratio2 = 0.05
  
readcol, 'previous/CH_Fisher_1995.txt', r, lgte
p30 = plot(r, lgte, ' o', over=p39, sym_size=rs_size, color=ct_a[*, i])
t30 = text(tpos2[0], tpos2[1]-ratio2*i, 'Fisher (1995)', /relative, target=p39, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)                                          
i++

readcol, 'previous/CH_Ko_1997.txt', r, lgte
p31 = plot(r, lgte, ' o', over=p39, sym_size=rs_size, color=ct_a[*, i])
t31 = text(tpos2[0], tpos2[1]-ratio2*i, 'Ko et al. (1997)', /relative, target=p39, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

readcol, 'previous/CH_David_1998.txt', r, lgte
p32 = plot(r, lgte, ' o', over=p39, sym_size=rs_size, color=ct_a[*, i])
t32 = text(tpos2[0], tpos2[1]-ratio2*i, 'David et al. (1998)', /relative, target=p39, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

tpos2 = [0.5, 0.23+ratio2*3.]

readcol, 'previous/CH_Wilhelm_2006.txt', r, lgte1, lgte2
p331 = plot(r, lgte1, ' o', over=p39, sym_size=rs_size, color=ct_a[*, i])
p332 = plot(r, lgte2, ' tu', over=p39, sym_size=rs_size, color=ct_a[*, i])
t32 = text(tpos2[0], tpos2[1]-ratio2*i, 'Wilhelm (2006)', /relative, target=p39, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

readcol, 'previous/CH_Landi_2008.txt', r, lgte
p34 = plot(r, lgte, ' o', over=p39, sym_size=rs_size, color=ct_a[*, i])
t34 = text(tpos2[0], tpos2[1]-ratio2*i, 'Landi (2008)', /relative, target=p39, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++

readcol, 'previous/CH_Reginald_2009.txt', r, lgte
p35 = plot(r, lgte, ' o', over=p39, sym_size=rs_size, color=ct_a[*, i])
t35 = text(tpos2[0], tpos2[1]-ratio2*i, 'Reginald et al. (2009)', /relative, target=p39, $
           color=ct_a[*, i], font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0)
i++


t36 = text(tpos2[0], tpos2[1]-ratio2*i, 'Region E', /relative, target=p39, $
           color='Black', font_size=rf_size, font_name='malgun gothic', $
           vertical_align=0.5, align=0, font_style=1)

dum1 = p19.title
t93 = text(0.25, 0.942, '(a) Region of Interest', $
            /normal, color='white', align=0.5, vertical_align=0., $
    font_size=12, font_style=1, font_name='malgun gothic')
;stop
w01.save, 'fig6.pdf', page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2, /bitmap


end