cd, '/data/home/chokh/tomography-1902'

if 1 then begin 
  restore, 'tdem_result.sav' ;; box_all
  tbin1 = 10.^(lgtmin+findgen(nlgt)*dlgt)
  tbin2 = 10.^(lgtmin+(findgen(nlgt)+1.)*dlgt)
  tbin = 0.5*(tbin1+tbin2)
;  nlgt = 12
;  tbin = tbin[0:nlgt-1]
  box_all = box_all[*, *, *, 0:nlgt-1]
  tbinn = rebin(reform(tbin, 1, 1, 1, nlgt), sz, sz, sz, nlgt)
  avg_te = total(tbinn*box_all, 4, /nan)/total(box_all, 4, /nan)
;  stop
  save, avg_te, tbin, sz, lgtmin, nlgt, dlgt, r_sun, $
        filename='avg_te.sav'
  stop
endif else restore, 'avg_te.sav'


; 2. Electron temperature

w01 = window(dim=[8d2, 8d2], background_color='black')

xp = findgen(128)-64
xxp = rebin(xp, 128, 128, 128)
yyp = rebin(reform(xp, 1, 128), 128, 128, 128)
zzp = rebin(reform(xp, 1, 1, 128), 128, 128, 128)
dist = sqrt(xxp^2.+yyp^2.+zzp^2.)
aa = dist*0.
aa[where(dist lt r_sun)] = 100
disk_in = volume(aa, /current, axis=0, $
                xr=[0, 128], yr=[0, 128], zr=[0, 128], $
                aspect_z=1, background_color='black', $
                pos=[0.2, 0.27, 0.8, 0.87])           
disk_in.hide = 1
disk_in.rotate, 75, /zaxis

restore, '/data/home/chokh/tomography/sdo_sinogram_0193.sav'
carr = where(t_obs ge -27.27*0.5 and t_obs lt 27.27*0.5)
surf_image = reverse(reform(sinogram[64, carr, 63.5-r_sun:63.5+r_sun]), 1)
surf_image = bytscl(surf_image, 0, 5d2)
aia_lct, rr, gg, bb, wave=193
ct193 = obj_new('IDLgrPalette', rr, gg, bb)
oImage = obj_new('IDLgrImage', surf_image, palette=ct193)
mesh_obj, 4, vert, poly, replicate(1, 101, 101)
vec = findgen(101)/100.
tex_coord = fltarr(2, 101, 101)
tex_coord[0, *, *] = rebin(vec, 101, 101)
tex_coord[1, *, *] = rebin(transpose(vec), 101, 101)
oSolar = obj_new('IDLgrPolygon', data=vert, polygons=poly, $
                  shading=1, color=[255, 255, 255], $
                  texture_coord=tex_coord, texture_map=oImage, $
                  /texture_interp)
surf = plot3d([63.5], [63.5], [63.5], over=disk_in, $
              sym_object=oSolar, $
              /data, sym_size=r_sun*0.615)
  
 

opa0 = findgen(256)*0.25
;opa0 = fltarr(256)+255
opa0[0] = 0

drange = [5.1d0, 7.0d0]
;rgbt = colortable(57, /rev)
rgbt = colortable(39)
v01 = volume(alog10(avg_te), axis=0, over=surf, render_qual=2, $
             rgb_table0=rgbt, opacity_table0=opa0, extents_trans=100, $
             min=drange[0], max=drange[1])
v02 = volume(alog10(avg_te), over=v01, render_qual=2, $
             rgb_table0=rgbt, opacity_table0=opa0, min=drange[0], max=drange[1])
v03 = volume(alog10(avg_te), over=v01, render_qual=2, $
             rgb_table0=rgbt, opacity_table0=opa0, min=drange[0], max=drange[1])

theta1 =  10.*!dtor
theta2 = 115.*!dtor
v01.clip_planes = -[tan(theta1), -1, 0, 64.*(1.-tan(theta1))]
v02.clip_planes = [tan(theta2), -1, 0, 64.*(1.-tan(theta2))]
v03.clip_planes = [0, 0, 1, -64]

pdum = plot(/test, /current, axis=0)
pdum.hide = 1
cb1 = colorbar(color='white', pos=[0.2, 0.1, 0.8, 0.12], /nor, $
               rgb_table=rgbt, range=drange, /border, $
               title='log(N!le!n) (cm!u-3!n)', $
               font_size=13, font_style=1, font_name='malgun gothic')

lcol = 'cyan'
lthick = 2

theta = [0:90:1]*!dtor
xdum = 64.5*cos(theta)*cos(theta1+!dpi)+63.5
ydum = 64.5*cos(theta)*sin(theta1+!dpi)+63.5
zdum = 64.5*sin(theta)+63.5
p10 = plot3d(xdum, ydum, zdum, color=lcol, thick=lthick, over=v01, clip=0)

xdum = 64.5*cos(theta)*cos(theta2)+63.5
ydum = 64.5*cos(theta)*sin(theta2)+63.5
zdum = 64.5*sin(theta)+63.5
p11 = plot3d(xdum, ydum, zdum, color=lcol, thick=lthick, over=v01, clip=0)

theta = [115:190]*!dtor
xdum = 64.5*cos(theta)+63.5
ydum = 64.5*sin(theta)+63.5
zdum = replicate(63.5, n_elements(theta))
p12 = plot3d(xdum, ydum, zdum, color=lcol, thick=lthick, over=v01, clip=0)

p13 = plot3d(63.5*[1, 1], 63.5*[1, 1], [r_sun, 64.5]+63.5, $
             color=lcol, thick=lthick, over=v01, clip=[0])

rdum = [r_sun+1, 64.5]
xdum = rdum*cos(theta[2])+63.5
ydum = rdum*sin(theta[2])+63.5
zdum = rdum*0.+64.5
p14 = plot3d(xdum, ydum, zdum, color=lcol, thick=lthick, over=v01, clip=0)

rdum = [r_sun, 64.5]
xdum = rdum*cos(theta[-3])+63.5
ydum = rdum*sin(theta[-3])+63.5
zdum = rdum*0.+64.5
p15 = plot3d(xdum, ydum, zdum, color=lcol, thick=lthick, over=v01, clip=0)

stop
w01.save, 'avg_te_fig1.png', resol=200, /bitmap


if 0 then begin
  cd, 'avg_te_slice'
  opa1 = fltarr(256)+200
  opa1[0] = 0
  w02 = window(dim=[8d2, 8d2], background_color='black')
  disk_in1 = volume(aa, /current, axis=0, $
                  xr=[0, 128], yr=[0, 128], zr=[0, 128], $
                  aspect_z=1, background_color='black', $
                  pos=[0.2, 0.27, 0.8, 0.87], $
                  rgb_table0=6, opacity_table0=opa1)
;  disk_in1.hide = 1
  disk_in1.rotate, 70, /zaxis
  v10 = volume(alog10(avg_te), axis=0, over=disk_in1, render_qual=2, $
               rgb_table0=rgbt, opacity_table0=opa0, extents_trans=100, $
               min=drange[0], max=drange[1])
  pdum = plot(/test, /current, axis=0)
  pdum.hide = 1
  cb2 = colorbar(color='white', pos=[0.2, 0.1, 0.8, 0.12], /nor, $
               rgb_table=rgbt, range=drange, /border, $
               title='log(N!le!n) (cm!u-3!n)', $
               font_size=13, font_style=1, font_name='malgun gothic')

  ct01 = obj_new('IDLgrPalette')
  ct01 -> loadct, 1
  buff = 2
  proc = 0
  for j=0, 128-1-proc do begin
;    j = 17
    v10.clip_planes = [-1, 0, 0, j]
    disk_in1.clip_planes = [-1, 0, 0, j]
    cross = fltarr(sz+buff*2, sz+buff*2)+255
    cross[buff:-1-buff, buff:-1-buff] $
          = bytscl(alog10(reform(avg_te[j+proc, *, *])>1.), drange[0], drange[1])

    oCrimg = obj_new('IDLgrImage', cross, palette=ct01)
    cross_pos = [[0, -64-buff, -64-buff], [0, 64+buff, -64-buff], $
                 [0, 64+buff, 64+buff], [0, -64-buff, 64+buff]]
    oCross = obj_new('IDLgrPolygon', cross_pos, color=[255, 255, 255], $
                     texture_coord=[[0, 0], [1, 0], [1, 1], [0, 1]], $
                     texture_map=oCrimg, alpha_channel=[0.7])
    plane10 = plot3d([j], [63.5], [63.5], over=disk_in1, $
                     sym_object=oCross, sym_size = 0.615, /data)                  
    t01 = text(j-1, 124, 5, /data, clip=0, $
               'Distance from the center of the Sun : '+$
               string((j+proc-63.5)/r_sun, f='(f5.2)')+' R$_\odot$', $
               font_size=13, font_style=1, font_name='malgun gothic', $
               baseline=[0, -1, 0], updir=[0, 0, 1], color='white')
;    stop               
    w02.save, string(j, f='(i03)')+'.png', resol=200
    t01.delete
    plane10.delete
  endfor
              
file = file_search('*.png')
ffmpeg, file, 20

cd, '..'
endif

end