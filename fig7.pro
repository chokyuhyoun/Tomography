path = '/hae/homedata/khcho/tomography-1906'
sav_path = path+'/k_corona_sav_5.5_6.5'
file_mkdir, sav_path

cd, path
restore, 'tdem_result_5.5_6.5.sav' ;; tot_e, avg_te
sz = (size(tot_e))[1]
xp = findgen(sz)-0.5*(sz-1.)
zp = findgen(sz)-0.5*(sz-1.)
xxp = rebin(xp, sz, sz)
zzp = rebin(transpose(zp), sz, sz)
dist = sqrt(xxp^2.+zzp^2.)
non_disk = where(dist gt r_sun, nreal)
coord = array_indices(dist, non_disk)

co_y = xp/r_sun

i_t = fltarr(128, 128)
i_r = fltarr(128, 128)
;stop
if 0 then begin
  cd, sav_path
  f = file_search('*.sav', count=dum1)
  if dum1 ne 0 then file_delete, f
  split_for, 0, nreal-1, commands=[$
    'density = tot_e[coord[0, i], *, coord[1, i]]', $
    'itir = k_corona_white(dist[non_disk[i]]/r_sun, density=density, $', $
    '                     co_z=co_y, i0=i0, /pol)', $
    'i_t[coord[0, i], coord[1, i]] = itir[0]', $
    'i_r[coord[0, i], coord[1, i]] = itir[1]'], $
    after_loop_command= $
       'save, i_t, i_r, i0, startp, endp, filename=string(i, f="(i04)")+".sav"', $
    varnames=['tot_e', 'coord', 'dist', 'non_disk', 'r_sun', 'co_y', $
              'i_t', 'i_r'], $
    ctvariable_name='i', nsplit=12
  
  sav_file = file_search('????.sav')
  i_t1 = i_t*0.
  i_r1 = i_r*0.
  for i=0, n_elements(sav_file)-1 do begin
    restore, sav_file[i]
     i_t1 = i_t1+i_t
     i_r1 = i_r1+i_r 
  endfor
  cd, '..'
  save, i_t1, i_r1, i0, filename='coronal_image.sav'
endif else restore, 'coronal_image.sav'

pb_image = (i_t1-i_r1)/i0*1d6
pb_image[where(dist/r_sun le 1.05)] = 0.

r = [1.1, 1.15, 1.2]
nr = n_elements(r)
th = findgen(361)*!dtor+0.5*!dpi
nth = n_elements(th)

m_r = r*r_sun
m_rr = rebin(m_r, nr, nth)
m_th = rebin(transpose(th), nr, nth)

m_x = m_rr*cos(m_th)+63.5
m_y = m_rr*sin(m_th)+63.5
mod_pro = interpolate(pb_image, m_x, m_y)


;; (a) K-Cor pB image
obs_file = '20190215_180517_kcor_l1.5_extavg.fts'
mreadfits, obs_file, h_obs, cor_obs
cor_obs = (cor_obs > 0.)*1d6
xp1 = (findgen(h_obs.naxis1)-h_obs.crpix1)/h_obs.r_sun
o_r = r*h_obs.r_sun
o_rr = rebin(o_r, nr, nth)
o_th = m_th
o_x = o_rr*cos(o_th)+h_obs.crpix1
o_y = o_rr*sin(o_th)+h_obs.crpix2
obs_pro = interpolate(cor_obs, o_x, o_y)

w01 = window(dim=[8d2, 8d2])
im01 = image_kh(cor_obs, xp1, xp1, /current, $
                pos=[100, 500, 350, 750], /dev, $
                min=0, max=5d-7*1d6, xthick=1.5, ythick=1.5, $
                xtitle='Solar X (R$_\odot$)', ytitle='Solar Y (R$_\odot$)', $
                xtickdir=1, ytickdir=1, xticklen=0.03, yticklen=0.03, $
                title='(a) K-Cor pB Image', font_size=12, $
                xr=minmax(co_y), yr=minmax(co_y))
p01 = plot(cos(findgen(361)*!dtor), sin(findgen(361)*!dtor), $
                over=im01, '-w', /data)

col = ['red', 'green', 'blue']
for i=0, nr-1 do begin
  p02 = plot(r[i]*cos(th), r[i]*sin(th), over=im01, '-', color=col[i])
endfor

cb01 = colorbar(target=im01, /relative, pos=[0, -0.30, 1, -0.25], $
                title='10$^{-6}$ B$_{Cor}$/B$_\odot$', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5)

;; (b) Synthetic coronal pB image
im02 = image_kh(pb_image, co_y, co_y, /current, $
                pos=[500, 500, 750, 750], /dev, $
                min=0, max=0.5, xthick=1.5, ythick=1.5, $
                xtitle='Solar X (R$_\odot$)', ytitle='Solar Y (R$_\odot$)', $
                xtickdir=1, ytickdir=1, xticklen=0.03, yticklen=0.03, $
                title='(b) Synthetic Coronal pB Image', font_size=12)

cb02 = colorbar(target=im02, /relative, pos=[0, -0.30, 1, -0.25], $
                title='10$^{-6}$ B$_{Cor}$/B$_\odot$', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5)


p03 = plot(cos(findgen(361)*!dtor), sin(findgen(361)*!dtor), $
           over=im02, '-w', /data)
for i=0, nr-1 do begin
  p04 = plot(r[i]*cos(th), r[i]*sin(th), over=im02, '-', color=col[i])
endfor

ang_txt = string(indgen(4)*90, f='(i0)')+'$\deg$'
for i=0, 3 do begin
  ang = (i+1.)*90*!dtor
  t01 = text(0.8*cos(ang), 0.8*sin(ang), ang_txt[i], $
              vertical_align=0.5, align=0.5, font_size=10, color='white', $
              target=im01, /data)
  t02 = text(0.8*cos(ang), 0.8*sin(ang), ang_txt[i], $
              vertical_align=0.5, align=0.5, font_size=10, color='white', $
              target=im02, /data)
endfor

p22 = plot(indgen(2), /nodata, pos=[100, 80, 750, 330], /dev, /current, $
           xr=[0, 360], yr=[0, 0.7], xmajor=5, xminor=2, yminor=1, $ 
           font_size=12, font_style=1, font_name='malgun gothic', $
           xtitle='Position Angle ($\deg$)', ytitle='Brightness ($10^{-6} B_{\odot})$', $
           xthick=1.5, ythick=1.5, $
           title='(c) Brightness Profiles')

p20 = objarr(nr)
p21 = objarr(nr)
for i=0, nr-1 do begin
  p20[i] = plot(findgen(361), mod_pro[i, *], over=p22, '-2', color=col[i])
  p21[i] = plot(findgen(361), obs_pro[i, *], over=p22, '--2', color=col[i])   
  t220 = text(0.85*p22.xr[1], (0.88-i*0.08)*p22.yr[1], /data, target=p22, $
              string(r[i], f='(f4.2)')+' R$_\odot$', $
              color=col[i], vertical_align=0.5, $
              font_size=12, font_style=1, font_name='malgun gothic')  
endfor
leg1 = string(r, f='(f4.2)')+' R$_\odot$'
p200 = plot([0.5, 0.55]*p22.xr[1], 0.88*[1, 1]*p22.yr[1], '-2', over=p22)
t200 = text(0.57*p22.xr[1], 0.88*p22.yr[1], 'Reconstructed', /data, $
            vertical_align=0.5, target=p22, $
            font_size=12, font_style=1, font_name='malgun gothic')
p210 = plot([0.5, 0.55]*p22.xr[1], 0.75*[1, 1]*p22.yr[1], '--2', over=p22)
t210 = text(0.57*p22.xr[1], 0.75*p22.yr[1], 'K-Cor', /data, $
            vertical_align=0.5, target=p22, $
            font_size=12, font_style=1, font_name='malgun gothic')
;stop

;w02.save, 'pa_profile1.png', resol=200
w01.save, 'fig7.pdf', page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2



end
