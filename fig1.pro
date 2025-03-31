path = '/data/home/chokh/tomography-1902'
cd, path

restore, 'tdem_result.sav' ;; box_all
restore, 'sdo_img_all.sav' ;; img

w01 = window(dim=[8d2, 8d2], background_color='black')

;; sphere -------------------------------------------------------------
;r_sun = 1.
range = 1.3

cast = volume(findgen(128, 128, 128), /current, /nodata, axis_style=2, $
              xr=range*[-1, 1], yr=range*[-1, 1], zr=range*[-1, 1], $
              aspect_z=1, aspect_ratio=1, pos=pos1, $
              xtickformat='(a1)', ytickformat='(a1)', ztickformat='(a1)', $
              xticklen=0, yticklen=0, zticklen=0, $
              xcolor='white', ycolor='white', zcolor='white')
pos1 = [0.35, 0.3, 0.85, 0.8]
cast.pos = pos1
mar=0.2
;ratio = (1.-2.*mar)/range
ratio = 67.*(pos1[2]-pos1[0])/range ;;empirical

;cast.hide=1
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
lct193 = [[rr], [gg], [bb]]
im91 = image(surf_img, findgen(im_sz[0])/im_sz[0]*360, $
              findgen(im_sz[1])/im_sz[1]*180.-90., $
              /current, pos=[0, 0, 1, 1], $
              rgb_table=lct193, min=0, max=4d2)
;stop

surf_data = w02.copywindow(border=0, resolution=300)
oImage = obj_new('IDLgrImage', shift(surf_data, 0, (size(surf_data))[2]*0.25, 0))

mesh_obj, 4, vert, poly, replicate(1, 101, 101)
vec = findgen(101)/100.
tex_coord = fltarr(2, 101, 101)
tex_coord[0, *, *] = rebin(vec, 101, 101)
tex_coord[1, *, *] = rebin(transpose(vec), 101, 101)
oSolar = obj_new('IDLgrPolygon', data=vert, polygons=poly, $
                  color=[255, 255, 255], $
                  texture_coord=tex_coord, texture_map=oImage, $
                  /texture_interp)
surf = plot3d([0], [0], [0], over=cast, $
              sym_object=oSolar, /data, sym_size=ratio)
cast.rotate, /reset
cast.rotate, -65, /x
cast.rotate, -50, /z
w02.close

;; arrow -------------------------------------------------------------
arrow_arr = fltarr(3, 5)
arrow_arr[0, *] = [0, 0.015, 0.015, 0.06, 0]
arrow_arr[2, *] = [-1.5, -1.5, 1.4, 1.4, 1.5]
mesh_obj, 6, vts, pol, arrow_arr, p1=100, /closed

vts = transpose(vts)
Oarrow1 = obj_new('IDLgrPolygon', transpose(vts), polygons=pol, $
                  color=[200, 200, 200])
arrow1 = plot3d([0], [0], [0], over=cast, clip=0, $
                sym_object=Oarrow1, /data, sym_size=ratio)
t_arrow = text(0, 0, 1.6, 'Rotation Axis', /data, target=cast, $
          color='white', align=0.5, vertical_align=0.5, clip=0, $
          font_size=18, font_style=1, font_name='malgun gothic', $
          updir=[0, 0, 1], baseline=[0, 1, 0])
radi = 0.2
ang = [!dpi : 2.3*!dpi : 0.01]+45.*!dtor
rx = radi*cos(ang)
ry = radi*sin(ang)
rot_arr1 = plot3d(rx, ry, replicate(range+0.05, n_elements(ang)), $
                 '-w2', over=cast, /data, clip=0)
endpx = rx[-1]+cos(ang[-1]+0.5*!dpi+[135., 0, 225.]*!dtor)*[0.1, 0, 0.1]
endpy = ry[-1]+sin(ang[-1]+0.5*!dpi+[135., 0, 225.]*!dtor)*[0.1, 0, 0.1]
rot_arr2 = plot3d(endpx, endpy, replicate(range+0.05, 3), '-w2', $
                  over=cast, /data, clip=0)

;; plane  -------------------------------------------------------------

lat_pos_a = 1.1
lat_pos_b = -0.2
rad_at_b = sqrt(1.-lat_pos_b^2)
plane_a = polygon([-range, range, range, -range], $
                  [-range, -range, range, range], $
                  replicate(lat_pos_a, 4), $
                  over=cast, /data, clip=0, fill_color='orange red', $
                  transp=40, linestyle='')
plane_b = polygon([-range, range, range, -range], $
                  [-range, -range, range, range], $
                  replicate(lat_pos_b, 4), $
                  over=cast, /data, clip=0, fill_color='blue', $
                  transp=40, linestyle='')

plane_b1 = polygon([-rad_at_b, rad_at_b, rad_at_b, -rad_at_b], $
                  [0, 0, range, range], $
                  replicate(lat_pos_b, 4)+0.001, $
                  over=cast, /data, clip=0, fill_color='white', $
                  pattern_orient=45, pattern_spac=7, pattern_thick=1, $
                  linestyle='')
t_plane_a = text(range-0.1, -range+0.4, lat_pos_a+0.01, target=cast, /data, $
                 'Plane 1', color='white', $
                 font_size=18, font_style=1, font_name='malgun gothic', $
                 align=0, vertical_align=0, $
                 updir=[-1, 0, 0], baseline=[0, 1, 0])
t_plane_b = text(range-0.1, -range+0.4, lat_pos_b+0.01, target=cast, /data, $
                 'Plane 2', color='white', $
                 font_size=18, font_style=1, font_name='malgun gothic', $
                 align=0, vertical_align=0, $
                 updir=[-1, 0, 0], baseline=[0, 1, 0])
;stop

;; observed image -----------------------------------------------------
buff = 1
obs_plane_pos = -2.5  ;; y dir
cross = fltarr(sz+buff*2, sz+buff*2)+255
cross[buff:-1-buff, buff:-1-buff] = bytscl(reform(img[*, *, 815, 3]), 0, 4d2)
oPalette = obj_new('IDLgrPalette', rr, gg, bb)
oCrimg = obj_new('IDLgrImage', cross, palette=oPalette)
cross_pos = [[-range, 0, -range], [range, 0, -range], $
             [range, 0, range], [-range, 0, range]]
oCross = obj_new('IDLgrPolygon', cross_pos, color=[255, 255, 255], $
                  texture_coord=[[0, 0], [1, 0], [1, 1], [0, 1]], $
                  texture_map=oCrimg, alpha_channel=1, $
                  linestyle=0)
obs_plane = plot3d([0], [obs_plane_pos], [0], over=cast, /data, clip=0, $
                   sym_object=oCross, sym_size=ratio)
t_obs_plane = text(0, obs_plane_pos, -range-0.1, $
                   'SDO/AIA 193 $\AA$ Image', $
                   target=cast, /data, clip=0, color='white', $
                   font_size=18, font_style=1, font_name='malgun gothic', $
                   align=0.5, vertical_align=0.5, $
                   updir=[0, 0, 1], baseline=[1, 0, 0]) 
p01 = plot3d([-range, range], obs_plane_pos*[1, 1], lat_pos_a*[1, 1], $
             over=cast, '-5r', clip=0)
p011 = plot3d(-range*[1, 1], [obs_plane_pos, -range], lat_pos_a*[1, 1], $
              over=cast, '--3r', clip=0)
p012 = plot3d(range*[1, 1], [obs_plane_pos, -range], lat_pos_a*[1, 1], $
              over=cast, '--3r', clip=0)

p02 = plot3d([-range, range], obs_plane_pos*[1, 1], lat_pos_b*[1, 1], $
             over=cast, '-5b', clip=0)
p021 = plot3d(-range*[1, 1], [obs_plane_pos, -range], lat_pos_b*[1, 1], $
             over=cast, '--3b', clip=0)
p022 = plot3d(range*[1, 1], [obs_plane_pos, -range], lat_pos_b*[1, 1], $
             over=cast, '--3b', clip=0)

;; axis

ar1 = [0.6, 0, 0, 1]
ar2 = [0.55, 0.05, 0, 1]
ar3 = [0.55, -0.05, 0, 1]
xrot = [0, 0, 0, 0, 0]
yrot = [0, 0, -90, 0, 0]
zrot = [0, 90, 0, -45, 45]

rp = [2.7, -0.8, -1.3]
axischa = ["!18x'", "y'", 'z', 'x', 'y']
t41 = objarr(n_elements(xrot))
for i=0, n_elements(xrot)-1 do begin
  t3d, /reset, rotate=[xrot[i], yrot[i], zrot[i]]
  ar11 = !p.t##ar1+rp
  ar22 = !p.t##ar2+rp
  ar33 = !p.t##ar3+rp
  p41 = plot3d([rp[0], ar11[0]], [rp[1], ar11[1]], [rp[2], ar11[2]], $
               over=cast, color='white', thick=2, clip=0)
  p42 = plot3d([ar11[0], ar22[0]], [ar11[1], ar22[1]], [ar11[2], ar22[2]], $
               over=cast, color='white', thick=2, clip=0)
  p43 = plot3d([ar11[0], ar33[0]], [ar11[1], ar33[1]], [ar11[2], ar33[2]], $
               over=cast, color='white', thick=2, clip=0)
  chpos = !p.t##(ar1+[0.15, 0, 0, 0])+rp
  t41[i] = text(chpos[0], chpos[1], chpos[2], axischa[i], over=cast, /data, $
                color='white', align=0.5, vertical_align=0.5, $
                clip=0, font_size=15, font_style=0, $
                updir=[cos(135*!dtor), sin(135*!dtor), 0], $
                baseline=[cos(45*!dtor), sin(45*!dtor), 0])
  t41[i].font_name='Hershey 8'
endfor
t41[2].updir=[0, 0, 1]
p44 = plot3d(0.3*cos(-findgen(45)*!dtor)+rp[0], $
             0.3*sin(-findgen(45)*!dtor)+rp[1], $
             replicate(rp[2], 45), $
             over=cast, /data, '-w2', clip=0)
t42 = text(0.5*cos(-22.5*!dtor)+rp[0], 0.5*sin(-22.5*!dtor)+rp[1], rp[2], $
           '$\theta$', color='white', over=cast, /data, $
           align=0.5, vertical_align=0.5, clip=0, font_size=15, font_style=1, $
          updir=t41[0].updir, baseline=t41[0].baseline)
             
 
w01.save, 'fig1.pdf', page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2, /bitmap

end