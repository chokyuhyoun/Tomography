path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path


w01 = window(dim=[8d2, 8d2])


; 1. Comparison between parallel beam and fan beam
restore, '3d_reconst_parallel.sav', /verbose
;stop
;nhead = rot_cube(head3d, [0, 0, 25])
nhead = head3d
sz = (size(head3d))[1]
xc = 0.5*(sz-1)
binsize = 1

obs_ypos = -1.5d8/(725.*0.6*32.)+xc   ;; in 32 binning pixel
obs_ypos = -sz

;; (a) Parallel Beam Illustration
pos01 = [0.1, 0.63, 0.33, 0.86]
mesh_obj, 4, vert, poly, replicate(r_sun, 201, 201)
oSolar = obj_new('IDLgrPolygon', data=vert, polygons=poly, $
                color=[255, 255, 255])
v01 = plot3d(findgen(sz), findgen(sz), findgen(sz), /nodata, $
              /current, axis=0, $
              xr=[0, 128], yr=[0, 128], zr=[0.01, 128], $
              aspect_ratio=1, aspect_z=1)
v01.pos = pos01
v01.rotate, -80, /x, /reset
v01.rotate, 50, /z
v02 = polygon(vert+xc, connectivity=poly, target=v01, /data, $
              fill_color='black', linestyle=6)

xray = 30.*(findgen(5)-2)+xc
zray = xray
nray = n_elements(xray)
xxray = rebin(xray, nray, nray)
zzray = rebin(transpose(zray), nray, nray)
dist1 = sqrt((xxray-xc)^2.+(zzray-xc)^2.)
y_end1 = 63.5-sqrt(r_sun^2.-dist1^2.)
y_end1[where(dist1 gt r_sun, /null)] = sz 
p01 = objarr(nray, nray)
for j=0, nray-1 do begin
  for i=0, nray-1 do begin
    p01[i, j] = plot3d(xxray[i, j]*[1, 1], [-xc, y_end1[i, j]], zzray[i, j]*[1, 1], $
      over=v01, /data, clip=0, color='cyan', thick=2)
  endfor
endfor
opa_table = findgen(256)*0.2
v03 = volume(nhead, over=v01, opacity_table0=opa_table)

;; (b) Parallel beam projection histogram
diff_par = box-head3d
real = where(finite(diff_par))
hist_par = histogram(diff_par[real], location=x1, binsize=binsize)
pos11 = [470, 470, 760, 760]
p11 = plot(x1, hist_par*1d-5, /histogram, pos=pos11, /dev, /current, $
          xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', $
          axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^5$)', $
          title='(b) Historam of Parallel Beam Reconst.', $
          font_size=12, font_style=1, font_name='malgun gothic', $
          xr=[-10, 10], yr=[0, 6], color='black')
p12 = plot([0, 0], [0, p11.yr[1]], over=p11, ':1')
t11 = text(6, 1.5, target=p11, /data, $
          '$\sigma$ = '+string(stddev(diff_par, /nan), f='(f4.1)'), $
          align=0.5, vertical_align=0.5, color=p11[0].color, $
          font_size=11, font_name='malgun gothic')

t01 = text(mean(p11.title.pos[[0, 2]])-0.5, mean(p11.title.pos[[1, 3]]), $
           '(a) Illustration of Parallel Beam', $
           align = 0.5, vertical_align=0.5, $
           font_size=12, font_style=1, font_name='malgun gothic')

;; (c) Fan Beam Projection illustration
pos21 = pos01-[0, 0.5, 0, 0.5]

cost = (-obs_ypos+xc)/sqrt((-obs_ypos+xc)^2.+dist1^2.)
pos_ang = atan(zzray-xc, xxray-xc)
dist2 = dist1*cost

y_init2 = obs_ypos
y_end2 = xc-r_sun*sin(acos(cost)+acos(dist2/r_sun))
y_end2[where(dist2 gt r_sun, /null)] = sz 
rho_init2 = (obs_ypos-y_init2)*dist1/(-obs_ypos+xc)
rho_end2 = (-obs_ypos+sz)*dist1/(-obs_ypos+xc)

rho_end21 = r_sun*cos(acos(cost)+acos(dist2/r_sun))
rho_end2[where(dist2 le r_sun, /null)] = rho_end21[where(dist2 le r_sun)] 
        
v21 = plot3d(findgen(sz), findgen(sz), findgen(sz), /nodata, $
              /current, axis=0, $
              xr=[0, 128], yr=[0, 128], zr=[0.01, 128], $
              aspect_ratio=1, aspect_z=1)
v21.pos = pos21
v21.rotate, -70, /x, /reset
v21.rotate, 50, /z              
v22 = polygon(vert+xc, connectivity=poly, target=v21, /data, $
              fill_color='black', linestyle=6)
p20 = objarr(nray, nray)        
for j=0, nray-1 do begin
  for i=0, nray-1 do begin
    p20[i, j] = plot3d([rho_init2[i, j], rho_end2[i, j]]*cos(pos_ang[i, j])+xc, $
                       [y_init2, y_end2[i, j]], $
                       [rho_init2[i, j], rho_end2[i, j]]*sin(pos_ang[i, j])+xc, $
                       over=v21, /data, clip=0, color='dark orange', thick=2)
  endfor
endfor
v23 = volume(nhead, over=v21, opacity_table0=opa_table)

;; (d) Fan beam projection histogram
restore, '3d_reconst_fan.sav', /verbose
diff_fan = box-head3d
real = where(finite(diff_fan))
hist_fan = histogram(diff_fan[real], location=x2, binsize=binsize)
pos31 = pos11-[0, 400, 0, 400]
p31 = plot(x2, hist_fan*1d-5, /histogram, pos=pos31, /dev, /current, $
          xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', $
          axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^5$)', $
          title='(d) Historam of Fan Beam Reconst.', $
          font_size=12, font_style=1, font_name='malgun gothic', $
          xr=[-10, 10], yr=[0, 6], color='black')
p32 = plot([0, 0], [0, p31.yr[1]], over=p31, ':1')
t31 = text(6, 1.5, target=p31, /data, $
          '$\sigma$ = '+string(stddev(diff_fan, /nan), f='(f4.1)'), $
          align=0.5, vertical_align=0.5, color=p31[0].color, $
          font_size=11, font_name='malgun gothic')

t21 = text(mean(p31.title.pos[[0, 2]])-0.5, mean(p31.title.pos[[1, 3]]), $
          '(c) Illustration of Fan Beam', $
          align = 0.5, vertical_align=0.5, $
          font_size=12, font_style=1, font_name='malgun gothic')

cd, '/data/home/chokh/tomography-1902'
w01.save, 'fig2_3dfanbeam.pdf', resol=300, /bitmap, $
  page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2


end