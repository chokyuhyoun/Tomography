path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path


w01 = window(dim=[8d2, 4d2])


; 1. Comparison between parallel beam and fan beam
restore, '3d_reconst_fan.sav', /verbose
;stop
;nhead = rot_cube(head3d, [0, 0, 25])
nhead = head3d
sz = (size(head3d))[1]
xc = 0.5*(sz-1)
binsize = 1

obs_ypos = -1.5d8/(725.*0.6*32.)+xc   ;; in 32 binning pixel
obs_ypos = -sz

;; (a) Parallel Beam Illustration
pos01 = [0.1, 0.16, 0.43, 0.81]
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

tilt_ang_arr = [0, 2, 4, 7]*1.
tilt_plane_col = ['red', 'dark orange', 'green', 'blue']

plane_x = [0, 0, sz, sz, 0]
plane_y = [0, sz, sz, 0, 0]
plane_z = [xc, xc, xc, xc, xc]
arrow_arr = fltarr(3, 5)
arrow_arr[0, *] = [0, 0.015, 0.015, 0.06, 0]*r_sun
arrow_arr[2, *] = [-1.5, -1.5, 1.4, 1.4, 1.5]*r_sun
mesh_obj, 6, vts, pol, arrow_arr, p1=300, /closed

v03 = objarr(n_elements(tilt_ang_arr))
v04 = objarr(n_elements(tilt_ang_arr))
for i=0, n_elements(tilt_ang_arr)-1 do begin
  t3d, /reset, rotate=[tilt_ang_arr[i], 0, 0]
  pos=!p.t[0:2, 0:2]##([[plane_x], [plane_y], [plane_z]]-xc)+xc
  v03[i] = polygon(pos[*, 0], pos[*, 1], pos[*, 2], $
                   over=v01, /data, fill_color=tilt_plane_col[i], $
                   fill_transp=50)
  arr_pos = !p.t[0:2, 0:2]##transpose(vts)+xc
  v04[i] = polygon(arr_pos[*, 0], arr_pos[*, 1], arr_pos[*, 2], $
                   connectivity=pol, over=v01, /data, clip=0, $
                   color=tilt_plane_col[i], transp=0, $
                   fill_color=tilt_plane_col[i], fill_transp=0)
endfor

opa_table = findgen(256)*0.2
v05 = volume(nhead, over=v01, opacity_table0=opa_table)

;; (b) histogram

f = file_search('3d_reconst_fan*.sav')
pos11 = [470, 70, 760, 360]
p10 = objarr(n_elements(f))
t10 = objarr(n_elements(f))
for i=0, n_elements(f)-1 do begin
  restore, f[i]
  if ~n_elements(tilt_ang) then tilt_ang = 0.
  diff_fan = box-head3d
  real = where(finite(diff_fan))
  hist_fan = histogram(diff_fan[real], location=x1, binsize=binsize)
  if i eq 0 then begin
    p10[i] = plot(x1, hist_fan*1d-5, /histogram, pos=pos11, /dev, /current, $
      xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', axis_style=2, $
      xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^5$)', $
      title='(b) Historam (Fan Beam + Tilt) ', $
      font_size=12, font_style=1, font_name='malgun gothic', $
      xr=[-10, 10], yr=[0, 6], color=tilt_plane_col[i])
  endif else begin
    p10[i] = plot(x1, hist_fan*1d-5, /histogram, over=p10[0], $
                  thick=1.5, color=tilt_plane_col[i])
  endelse
  t10[i] = text(1, 5.2-0.5*i, target=p10[0], /data, $
                '$\sigma (\gamma$ ='+string(tilt_ang, f='(i0)')+'$\deg$) = '$
                +string(stddev(diff_fan, /nan), f='(f4.1)'), $
                align=0., vertical_align=0., color=p10[i].color, $
                font_size=11, font_name='malgun gothic')  
                
endfor

t01 = text(mean(p10[0].title.pos[[0, 2]])-0.5, mean(p10[0].title.pos[[1, 3]]), $
  '(a) Illustration of Axial Tilt of the Sun', $
  align = 0.5, vertical_align=0.5, $
  font_size=12, font_style=1, font_name='malgun gothic')
  
cd, '/data/home/chokh/tomography-1902'
w01.save, 'fig2_3dfanbeam_incli.pdf', resol=300, /bitmap, $
  page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2


end