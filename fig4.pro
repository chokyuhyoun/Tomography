path='/hae/homedata/khcho/tomography-1906'
cd, path

restore, 'sdo_img_all.sav' ; img
restore, 'sdo_indices.sav' ; oindex

t0 = julday(6, 15, 2019, 0, 0, 0)
t_obs = anytim2jul(oindex[*, 0].t_obs)
delt = t_obs-t0
sz=128

index = oindex[3, 0] ;; 193 A  1st  index
xc = index.crpix1*sz/index.naxis1
yc = index.crpix2*sz/index.naxis2
r_sun = index.r_sun*sz/index.naxis1

rot_co=[14.713, -2.396, -1.787]  ;; differential rotation coeff.

;; tomography
lat=asin((findgen(sz)-yc)/r_sun)
lat[where(~finite(lat) and (findgen(sz) lt xc))]=-0.5*!dpi
lat[where(~finite(lat) and (findgen(sz) gt xc))]=0.5*!dpi
ang_vel=rot_co[0]+rot_co[1]*sin(lat)^2.+rot_co[2]*sin(lat)^4.  ;  deg/day

if 1 then begin
  ti = systime(/sec)
  box=fltarr(sz, sz, sz)
  for j=0l, sz-1 do begin
    print, j
    theta=delt*ang_vel[j]
    box[*, *, j]=solar_tomography(reform(img[*, j, *, 3]), $
                                  theta, r_sun=r_sun, lat=lat[j])
  endfor
  print, (systime(/sec)-ti)
  save, box, filename='box_193.sav'
endif else restore, 'box_193.sav'

theta = -(rot_co[0]+rot_co[1]*sin(0)^2.+rot_co[2]*sin(0)^4.)*delt ;; deg
caldat, t_obs, mon, day, yr, hr, min, sec
mon = string(mon, f='(i2)')
day = string(day, f='(i2)')
yr = string(yr, f='(i4)')
hr = string(hr, f='(i02)')
min = string(min, f='(i02)')
;stop

xp=findgen(sz)-xc
xxp=rebin(xp, sz, sz, sz)
yyp=rebin(reform(xp, 1, sz), sz, sz, sz)
zzp=rebin(reform(xp, 1, 1, sz), sz, sz, sz)
dist=sqrt(xxp^2.+yyp^2.+zzp^2.)
inside=where(dist lt r_sun)
box[inside]=0
xxp1=rebin(xp, sz, sz)
yyp1=rebin(reform(xp, 1, sz), sz, sz)
dist1=sqrt(xxp1^2.+yyp^2.)
occ=rebin(reform(dist1, sz, 1, sz), sz, sz, sz)
filter = xxp*0+1
filter[where(occ lt r_sun and yyp gt 0)] = 0
ref_ang = 0. ; in deg
ref_t = arr_eq(theta, ref_ang)

aia_lct, rr, gg, bb, wave=193, /load
ct = [[rr], [gg], [bb]]
min = 0
max = 5d2
xp = (findgen(sz)-0.5*(sz-1))/r_sun

w01 = window(dim=[8d2, 8d2])
obs_img = reform(img[*, *, ref_t, 3])
im01 = image_kh(obs_img, xp, xp, /current, /dev, $
               pos=[100, 510, 350, 760], rgb_table=ct, $
               min=min, max=max, xthick=1.5, ythick=1.5, font_size=12, $
               xtitle='Solar X (R$_\odot$)', ytitle='Solar Y (R$_\odot$)', $
               title='(a) Observed Image', $
               xtickdir=1, ytickdir=1, xticklen=0.03, yticklen=0.03)
xx = cos(findgen(361)*!dtor)
yy = sin(findgen(361)*!dtor)

r1 = 0.2
xc1 = -0.95
yc1 = 0.15

r2 = 0.15
xc2 = 0.95
yc2 = -0.25

r3 = 0.2
xc3 = 1.05
yc3 = 0.25

p011 = plot(xx*r1+xc1, yy*r1+yc1, '-2', color='gray', over=im01, /data)
p012 = plot(xx*r2+xc2, yy*r2+yc2, '-2', color='gray', over=im01, /data)
p013 = plot(xx*r3+xc3, yy*r3+yc3, '-2', color='gray', over=im01, /data)

xc4 = 0.
yc4 = -sin(20.*!dtor)
xc5 = 0.
yc5 = -1.
t011 = text(xc1+(r1+0.08), yc1, 'A', target=im01, /data, $
            align=0.5, vertical_align=0.5, color='white', $
            font_style=1, font_size=9, font_name='malgun gothic')
t012 = text(xc2-(r2+0.08), yc2, 'B', target=im01, /data, $
            align=0.5, vertical_align=0.5, color='white', $
            font_style=1, font_size=9, font_name='malgun gothic')
t013 = text(xc3-(r3+0.08), yc3, 'C', target=im01, /data, $
            align=0.5, vertical_align=0.5, color='white', $
            font_style=1, font_size=9, font_name='malgun gothic')
t014 = text(xc4, yc4, 'D', target=im01, /data, $
            align=0.5, vertical_align=0.5, color='white', $
            font_style=1, font_size=9, font_name='malgun gothic')
t015 = text(xc5, yc5, 'E', target=im01, /data, $
            align=0.5, vertical_align=0.5, color='white', $
            font_style=1, font_size=9, font_name='malgun gothic')


cb01 = colorbar(target=im01, /relative, pos=[0, -0.28, 1, -0.23], $
               title='DN', orient=0, textpos=0, $
               font_style=1, font_size=11, font_name='malgun gothic', $
               /border, ticklen=0.5, minor=1, thick=1.5)
               
;box = rot_cube(box, [0, 0, 180])
if ref_ang ne 0 then box1 = rot_cube(box, [0, 0, ref_ang]) else box1 = box

box1[where(occ lt r_sun and yyp gt 0)]=0.
recon_img = total(box1, 2, /nan)
im02 = image_kh(recon_img, xp, xp, /current, /dev, $
               pos=[500, 510, 750, 760], rgb_table=ct, $
               min=min, max=max, xthick=1.5, ythick=1.5, font_size=12, $
               xtitle='Solar X (R$_\odot$)', ytitle='Solar Y (R$_\odot$)', $
               title='(b) Synthetic Image', $
               xtickdir=1, ytickdir=1, xticklen=0.03, yticklen=0.03)

cb02 = colorbar(target=im02, /relative, pos=[0, -0.28, 1, -0.22], $
               title='DN', orient=0, textpos=0, $
               font_style=1, font_size=11, font_name='malgun gothic', $
               /border, ticklen=0.5, minor=1, thick=1.5)

p021 = plot(xx*r1+xc1, yy*r1+yc1, '-2', color='gray', over=im02, /data)
p022 = plot(xx*r2+xc2, yy*r2+yc2, '-2', color='gray', over=im02, /data)
p023 = plot(xx*r3+xc3, yy*r3+yc3, '-2', color='gray', over=im02, /data)

dif_img = recon_img-obs_img
im03 = image_kh(dif_img, xp, xp, /current, /dev, $
                pos=[100, 110, 350, 360], rgb_table=33, $
                min=-100, max=100, xthick=1.5, ythick=1.5, font_size=12, $
                xtitle='Solar X (R$_\odot$)', ytitle='Solar Y (R$_\odot$)', $
                title='(c) Intensity Difference Image', $
                xtickdir=1, ytickdir=1, xticklen=0.03, yticklen=0.03)


p031 = plot(xx*r1+xc1, yy*r1+yc1, '-2', color='gray', over=im03, /data)
p032 = plot(xx*r2+xc2, yy*r2+yc2, '-2', color='gray', over=im03, /data)
p033 = plot(xx*r3+xc3, yy*r3+yc3, '-2', color='gray', over=im03, /data)

cb03 = colorbar(target=im03, /relative, pos=[0, -0.28, 1, -0.23], $
                title='$\Delta$ DN', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5)

hist_range = 50.*[-1, 1]
hist=histogram(dif_img, location=x1, min=hist_range[0], max=hist_range[1], $
               binsize=1)
p01=plot(x1, hist, pos=[480, 90, 750, 360], /dev, /current, $
  xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', ymajor=5, $
  axis_style=2, xtitle='$\Delta$ DN', ytitle='# of pixels', $
  title='(d) Intensity Difference Histogram', $
  font_size=12, font_style=1, font_name='malgun gothic', $
  /histogram, xr=hist_range, yr=[0, 2000])
p02=plot([0, 0], [0, p01.yr[1]], over=p01, ':1')
;fitres = mpfitpeak(x1+0.5, hist, arg, nterms=3, /lorentzian)
;p03 = plot(x1+0.5, fitres, '--r', over=p01)               
t03 = text(0.75, 0.5, /relative, target=p01,  $
          '$\sigma$ = '+string(stddev(dif_img, /nan), f='(f4.1)'), $
          align=0.5, vertical_align=0.5, $
          font_size=12, font_name='malgun gothic') 
;stop
w01.save, 'fig4.pdf', resol=300, $
            page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2  


if 0 then begin
res_path=path+'/img1'
file_mkdir, res_path
cd, res_path
num=0
w02=window(dim=[8d2, 4d2], background_color='black')
z_ang = 45.*sin(findgen(181)*!dtor*2.) 
i = 0
num=num+1
box1=rot_cube(box, [0, 0, z_ang[i]])
box1[where(occ lt r_sun and yyp gt 0)]=0.
t_arr = arr_eq(theta, z_ang[i])
im1 = image(reform(img[*, *, t_arr, 3]), pos=[0, 0, 400, 400], /dev, $
            /current, min=min, max=max, rgb_table=ct)
im2  =image(total(box1, 2, /nan), pos=[400, 0, 800, 400], /dev, $
            /current, min=min, max=max, rgb_table=ct)
t1=text(20, 20, 'SDO/AIA 193 $\AA$ '+strmid(oindex[t_arr, 3].t_obs, 0, 19), $
        /dev, color='white', font_size=15)
t2=text(430, 20, 'Rotation Angle : (0$\deg$, 0$\deg$, '+ $
        string(z_ang[i], f='(i0)')+'$\deg$)', $
        /dev, color='white', font_size=15)
;stop
w02.save, string(num, f='(i03)')+'.png', resol=200
if 1 then begin
for i=0, n_elements(z_ang)-1 do begin
  num=num+1
  box1=rot_cube(box, [0, 0, z_ang[i]])
  box1[where(occ lt r_sun and yyp gt 0)]=0.
  t_arr = arr_eq(theta, z_ang[i])

  setdata, im1, reform(img[*, *, t_arr, 3])
  setdata, im2, total(box1, 2, /nan)
  t1.string = 'SDO/AIA 193 $\AA$ '+strmid(oindex[t_arr, 3].t_obs, 0, 19)
  t2.string = 'Rotation Angle : (0$\deg$, 0$\deg$, '+ $
        string(z_ang[i], f='(i0)')+'$\deg$)'
  w02.save, string(num, f='(i03)')+'.png', resol=200
endfor
endif

stop
num = n_elements(z_ang)
for i=0, n_elements(z_ang)-1 do begin
  num=num+1
  box1=rot_cube(box, [z_ang[i], 0, 0])
  box1[where(occ lt r_sun and yyp gt 0)]=0.
  setdata, im2, total(box1, 2, /nan)
  t2.string = 'Rotation Angle : ('+string(z_ang[i], f='(i0)')+ $
               '$\deg$, 0$\deg$, 0$\deg$)'
  w02.save, string(num, f='(i03)')+'.png', resol=200
endfor
f = file_search('*.png')
ffmpeg, f, 20, filename='fig2_animation.mp4'
endif
end
