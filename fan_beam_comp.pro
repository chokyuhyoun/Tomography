path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path


w01 = window(dim=[8d2, 4d2])

minv = 0.
maxv = 255.
dminv = -20.
dmaxv = 20.

; 1. Comparison between parallel beam and fan beam
restore, 'sinogram1.sav', /verbose
;stop
rotang = 30.*!dtor
mx = [cos(rotang), -sin(rotang)]
my = [sin(rotang), cos(rotang)]
pos01 = [70, 50, 370, 350]

sz=(size(sinogram))[1]
xc=(sz-1)*0.5
yc=(sz-1)*0.5
xp=findgen(sz)-xc
yp=findgen(sz)-yc
xxp=rebin(xp, sz, sz)
yyp=rebin(transpose(yp), sz, sz)
data[where(sqrt(xxp^2.+yyp^2.) lt r_sun)] = 255.

hsize = 140
xcen = mean(pos01[[0, 2]])
ytop = pos01[3]
pos31 = [xcen-hsize, ytop-2.*hsize, xcen+hsize, ytop]
im31 = image(data, xp-0.5, yp-0.5, pos=pos31, /current, /dev, min=minv, max=maxv, $
  xr=minmax(xp), yr=minmax(yp), axis=0)
dashx1 = [(findgen(15)+1)*8]-xc
dashx2 = [(findgen(15)+1)*8]-yc
dashy1 = replicate(-200, n_elements(dashx1))
dashy2 = -sqrt(r_sun^2.-dashx1^2.)
dashy2[where(~finite(dashy2))] = 200.
dashx11 = mx##[[dashx1], [dashy1]]
dashy11 = my##[[dashx1], [dashy1]]
dashx22 = mx##[[dashx2], [dashy2]]
dashy22 = my##[[dashx2], [dashy2]]

p31 = objarr(n_elements(dashx1))
for i=0, n_elements(dashx1)-1 do begin
  p31[i] = plot([dashx11[i], dashx22[i]], [dashy11[i], dashy22[i]], $
                 over=im31, /data, $
                 '1', color='dark orange', clip=1)
endfor     


xp_pos = dashx1 ; x position at y = 63.5
slope = 1.5d8/(dashx1*725.*19.2*10.)
dashy3 = replicate(-200, n_elements(dashx1))
dashx3 = dashy3/slope+xp_pos
dashy4 = replicate(200, n_elements(dashx1))
dashx4 = dashy4/slope+xp_pos

ll = abs(-xp_pos*slope)/sqrt(slope^2.+1.) ;; distance from solar center to los line
sol_surfy = -r_sun*sin(acos(ll/r_sun)+0.5*!dpi-abs(atan(slope)))
sol_surfx = -r_sun*cos(acos(ll/r_sun)+0.5*!dpi+atan(slope))
dashx4[where(ll lt r_sun)] = sol_surfx[where(ll lt r_sun)]
dashy4[where(ll lt r_sun)] = sol_surfy[where(ll lt r_sun)]

;y_fan = yyp
;x_fan = (yyp-xc)/slope1+xp_pos       
;
;ll = abs((xc-xp)*slope)/sqrt(slope^2.+1.) ;; distance from solar center to los line
;sol_surf = -r_sun*sin(acos(ll/r_sun)+0.5*!dpi-abs(atan(slope)))+63.5
;sol_surf[where(~finite(sol_surf))] = 200.


dashx33 = mx##[[dashx3], [dashy3]]
dashy33 = my##[[dashx3], [dashy3]]
dashx44 = mx##[[dashx4], [dashy4]]
dashy44 = my##[[dashx4], [dashy4]]

p32 = objarr(n_elements(dashx1))
for i=0, n_elements(dashx1)-1 do begin
  p32[i] = plot([dashx33[i], dashx44[i]], [dashy33[i], dashy44[i]], $
                 over=im31, /data, $
                 '1', color='cyan', clip=1) 
endfor
arr = [75, 90, 100]
arx = arr*cos(1.5*!dpi+rotang)
ary = arr*sin(1.5*!dpi+rotang)
a31 = arrow(arx[0:1], ary[0:1], target=im31, clip=0, /data, $
            thick=1.5)
t31 = text(arx[2], ary[2],'Observer', /data, target=im31, clip=0, $
          font_size=12, font_style=0, font_name='malgun gothic', $
          align=0.5, orient=rotang*!radeg)

t05 = text(0, 0, 'Solar !cInterior', /data, $
  font_size=12, font_style=0, font_name='malgun gothic', $
  align=0.5, vertical_align=0.5, target=im31)

t06 = text(mean(pos31[[0, 2]]), pos31[3]+12, /dev, $
    '(a) Illustration of Fan Beam Reconstruction', $
    font_size=12, font_style=1, font_name='malgun gothic', $
    align=0.5, vertical_align=0)
t07 = text(-50, -74, 'Parallel Beam', target=im31, /data, clip=0, $
           color=p31[0].color, $
           font_size=12, font_style=0, font_name='malgun gothic')
t08 = text(-50, -84, 'Fan Beam (x10)', target=im31, /data, clip=0, $
           color=p32[0].color, $
           font_size=12, font_style=0, font_name='malgun gothic')

restore, 'solar_tomo_result.sav', /verbose
diff=res-data
real=where(finite(res))
hist=histogram(diff[real], location=x1, binsize=binsize)
pos21 = pos01+[400, 0, 400, 0]
p21=plot(x1, hist*1d-3, pos=pos21, /dev, /current, $
        xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', ymajor=5, $
        axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^3$)', $
        title='(b) Historam of Difference Image', $
        font_size=12, font_style=1, font_name='malgun gothic', $
        /histogram, xr=[-20, 20], yr=[0, 5], color=p31[0].color)
;stopo
p22=plot([0, 0], [0, p21.yr[1]], over=p21, ':1')
t21 = text(10, 1.5, target=p21, /data, $
          '$\sigma$ = '+string(stddev(diff, /nan), f='(f4.1)'), $
          align=0.5, vertical_align=0.5, color=p21[0].color, $
          font_size=11, font_name='malgun gothic')

restore, 'solar_tomo_fanbeam_result.sav', /verbose
diff=res-data
real=where(finite(res))
hist=histogram(diff[real], location=x1, binsize=binsize)

p23=plot(x1, hist*1d-3, over=p21, $
        thick=1.5, /histogram, color=p32[0].color)
t23 = text(10, 1.25, target=p23, /data, $
        '$\sigma$ = '+string(stddev(diff, /nan), f='(f4.1)'), $
        align=0.5, vertical_align=0.5, color=p23[0].color, $
        font_size=11, font_name='malgun gothic')


cd, '/data/home/chokh/tomography-1902'
w01.save, 'fig2_fanbeam.pdf', resol=300, /bitmap, $
  page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2
end