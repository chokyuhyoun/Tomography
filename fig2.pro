path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path


w01 = window(dim=[8d2, 8d2])

minv = 0.
maxv = 255.
dminv = -20.
dmaxv = 20.

; 1. Classical FBP

restore, 'tomo_result.sav', /verbose

diff=res-data
real=where(finite(res))
binsize = 1
hist=histogram(diff[real], location=x1, binsize=binsize)
pos01 = [70, 450, 370, 750]
p01=plot(x1, hist*1d-3, pos=pos01, /dev, /current, $
  xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', ymajor=5, $
  axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^3$)', $
  title='(a) Typical Sinogram', $
  font_size=12, font_style=1, font_name='malgun gothic', $
  /histogram, xr=[dminv, dmaxv], yr=[0, 5])
p02=plot([0, 0], [0, p01.yr[1]], over=p01, ':1')
;fitres = mpfitpeak(x1+0.5*binsize, hist*1d-3, arg, nterms=3, /lorentzian)
;p03 = plot(x1+0.5*binsize, fitres, '--r', over=p01)

pos021 = [20, 165, 115, 280]+pos01[[0, 1, 0, 1]]
pos022 = pos021+[0, -135, 0, -135]
pos023 = pos021+[160, 0, 160, 0]


im01=image(data, /current, title='Original', min=minv, max=maxv, $
  pos=pos021, /dev, font_size=10)
im02=image(res, /current, title='Reconstructed', min=minv, max=maxv, $
  pos=pos022, /dev, font_size=10)
im03=image(diff, /current, title='Difference', min=dminv , max=dmaxv, $
  pos=pos023, /dev, rgb_table=33, font_size=10)
c01=colorbar(target=im01, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=[minv, 0.5*(minv+maxv), maxv], textpos=0, minor=0)
c02=colorbar(target=im02, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c01.tickv, textpos=0, minor=0)
c03=colorbar(target=im03, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=[dminv, 0., dmaxv], textpos=0, minor=0)

t03 = text(mean(pos023[[0, 2]]), mean(pos022[[1, 3]]), /dev, $
           '$\sigma$ = '+string(stddev(diff, /nan), f='(f4.1)'), $
           align=0.5, vertical_align=0.5, $
           font_size=11, font_name='malgun gothic')     
;stop

; 2. Environment of solar observation
restore, 'sinogram1.sav', /verbose
;stop
rotang = 30.*!dtor
mx = [cos(rotang), -sin(rotang)]
my = [sin(rotang), cos(rotang)]

sz=size(sinogram)
xc=(sz[1]-1)*0.5
yc=(sz[1]-1)*0.5
xp=findgen(sz[1])-xc
yp=findgen(sz[1])-yc
xxp=rebin(xp, sz[1], sz[1])
yyp=rebin(transpose(yp), sz[1], sz[1])
data[where(sqrt(xxp^2.+yyp^2.) lt r_sun)] = 255.

hsize = 140
xcen = mean(pos01[[0, 2]])+400
ytop = pos01[3]
pos31 = [xcen-hsize, ytop-2*hsize, xcen+hsize, ytop]
im31 = image(data, xp, yp, pos=pos31, /current, /dev, min=minv, max=maxv, $
             xr=minmax(xp), yr=minmax(yp), axis=0)
dashx1 = [-r_sun*[1, 1], r_sun*[1, 1]]
dashy1 = [-200, 200, -200, 200]
dashx = mx##[[dashx1], [dashy1]]
dashy = my##[[dashx1], [dashy1]]
p31 = plot(dashx[0:1], dashy[0:1], over=im31, $
  '--2', color='yellow', clip=1)
p32 = plot(dashx[2:3], dashy[2:3], over=im31, $
  '--2', color='yellow', clip=1)

hatch = polygon([dashx[1], mean(dashx[0:1]), mean(dashx[2:3]), dashx[3]], $
                [dashy[1], mean(dashy[0:1]), mean(dashy[2:3]), dashy[3]], $
                target=im31, /data, clip=1, linestyle=6, $
                pattern_orient=45, fill_color='red', pattern_spacing=10)

cirang = [0:360:1]*!dtor
cirx = r_sun*cos(cirang)
ciry = r_sun*sin(cirang)

crosshatch1 = polygon(cirx, ciry, target=im31, /data, $
                      linestyle=6, fill_color='red', $
                      pattern_orient=45, pattern_spacing=10)
crosshatch2 = polygon(cirx, ciry, target=im31, /data, $
                      linestyle=6, fill_color='red', $
                      pattern_orient=-45, pattern_spacing=10)
             

arr = [80, 95, 105]
arx = arr*cos(1.5*!dpi+rotang)
ary = arr*sin(1.5*!dpi+rotang)
a31 = arrow(arx[0:1], ary[0:1], target=im31, clip=0, /data, $
  thick=1.5)
t31 = text(arx[2], ary[2],'Observer', /data, target=im31, clip=0, $
  font_size=12, font_style=0, font_name='malgun gothic', $
  align=0.5, orient=rotang*!radeg)

regr = 37.
regrotang = [1., 1.5, 0, 0.5]*!dpi
regx = regr*cos(regrotang+rotang)
regy = regr*sin(regrotang+rotang)

t32 = text(regx[0], regy[0], '(1)', /data, color='yellow',$
  font_size=12, font_style=0, font_name='malgun gothic', $
  align=0.5, vertical_align=0.5, target=im31)

t33 = text(regx[1], regy[1], '(2)', /data, color='yellow',$
  font_size=12, font_style=0, font_name='malgun gothic', $
  align=0.5, target=im31)

t34 = text(regx[2], regy[2], '(3)', /data, color='yellow',$
  font_size=12, font_style=0, font_name='malgun gothic', $
  align=0.5, vertical_align=0.5, target=im31)

t35 = text(regx[3], regy[3], '(4)', /data, color='yellow',$
    font_size=12, font_style=0, font_name='malgun gothic', $
    align=0.5, target=im31)

t05 = text(0, 0, 'Solar !cInterior', /data, $
  font_size=12, font_style=0, font_name='malgun gothic', $
  align=0.5, vertical_align=0.5, target=im31)

xcc = -25.
ycc = -95.
axisrotang = [0, 0.5, 0, 0.5]*!dpi+rotang*[0, 0, 1, 1]
axisr = 20
axisarr = objarr(4)
axistext = objarr(4)
axischa = ['!18x', 'y', "x'", "y'!3"]
for i=0, 3 do begin
  axisarr[i] = arrow([0, axisr*cos(axisrotang[i])]+xcc, $
                      [0, axisr*sin(axisrotang[i])]+ycc, $
                      target=im31, /data, clip=0, head_size=0.5)
  axistext[i] = text((axisr+5.)*cos(axisrotang[i])+xcc, $
                  (axisr+5.)*sin(axisrotang[i])+ycc, $
                  axischa[i], target=im31, clip=0, /data, $
                  align=0.5, vertical_align=0.5, $
                  font_size=9, font_style=0)
  axistext[i].font_name='Hershey 8'                   
endfor
angplot = plot(10.*cos([0:rotang:0.1])+xcc, 10.*sin([0:rotang:0.1])+ycc, $
               over=im31, /data, clip=0)
angtext = text(xcc+10, ycc-5, '$\theta$', target=im31, /data, clip=0, $
               align=0.5, vertical_align=0.5, $
               font_size=11, font_style=0, font_name='Hershey 18')


t06 = text(mean(pos31[[0, 2]]), pos31[3]+12, /dev, $
           '(b) Illustration of the Lat. Plane', $
           font_size=12, font_style=1, font_name='malgun gothic', $
           align=0.5, vertical_align=0)


; 3. Classical FBP with solar observation

restore, 'tomo_result1.sav', /verbose

diff=res-data
real=where(finite(res))
hist=histogram(diff[real], location=x1, binsize=binsize)
pos11 = pos01+[0, -400, 0, -400]
p11=plot(x1, hist*1d-3, pos=pos11, /dev, /current, $
  xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', ymajor=5, $
  axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^3$)', $
  title='(c) Sinogram from Solar EUV Obs.', $
  font_size=12, font_style=1, font_name='malgun gothic', $
  /histogram, xr=p01.xr, yr=p01.yr)
p12=plot([0, 0], [0, p11.yr[1]], over=p11, ':1')
;fitres = mpfitpeak(x1+0.5*binsize, hist*1d-3, arg, nterms=3, /lorentzian)
;p13 = plot(x1+0.5*binsize, fitres, '--r', over=p11)


pos121 = [20, 165, 115, 280]+pos11[[0, 1, 0, 1]]
pos122 = pos121+[0, -135, 0, -135]
pos123 = pos121+[160, 0, 160, 0]

im11=image(data, /current, title='Original', min=minv, max=maxv, $
  pos=pos121, /dev, font_size=10)
im12=image(res, /current, title='Reconstructed', min=minv, max=maxv, $
  pos=pos122, /dev, font_size=10)
im13=image(diff, /current, title='Difference', min=dminv, max=dmaxv, $
  pos=pos123, /dev, rgb_table=33, font_size=10)
c11=colorbar(target=im11, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c01.tickv, textpos=0, minor=0)
c12=colorbar(target=im12, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c01.tickv, textpos=0, minor=0)
c13=colorbar(target=im13, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c03.tickv, textpos=0, minor=0)

t13 = text(mean(pos123[[0, 2]]), mean(pos122[[1, 3]]), /dev, $
    '$\sigma$ = '+string(stddev(diff, /nan), f='(f4.1)'), $
    align=0.5, vertical_align=0.5, $
    font_size=11, font_name='malgun gothic')

; 4. Modified FBP with solar observation

restore, 'solar_tomo_result.sav', /verbose

diff=res-data
real=where(finite(res))
hist=histogram(diff[real], location=x1, binsize=binsize)
pos21 = pos01+[400, -400, 400, -400]
p21=plot(x1, hist*1d-3, pos=pos21, /dev, /current, $
  xthick=1.5, ythick=1.5, thick=1.5, ytickformat='(i0)', ymajor=5, $
  axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixels ($\times 10^3$)', $
  title='(d) Modified Sinogram', $
  font_size=12, font_style=1, font_name='malgun gothic', $
  /histogram, xr=p01.xr, yr=p01.yr)
p22=plot([0, 0], [0, p21.yr[1]], over=p21, ':1')
;fitres = mpfitpeak(x1+0.5*binsize, hist*1d-3, arg, nterms=3, /lorentzian)
;p23 = plot(x1+0.5*binsize, fitres, '--r', over=p21)


pos221 = [20, 165, 115, 280]+pos21[[0, 1, 0, 1]]
pos222 = pos221+[0, -135, 0, -135]
pos223 = pos221+[160, 0, 160, 0]

im21=image(data, /current, title='Original', min=minv, max=maxv, $
  pos=pos221, /dev, font_size=10)
im22=image(res, /current, title='Reconstructed', min=minv, max=maxv, $
  pos=pos222, /dev, font_size=10)
im23=image(diff, /current, title='Difference', min=dminv, max=dmaxv, $
  pos=pos223, /dev, rgb_table=33, font_size=10)
c21=colorbar(target=im21, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c01.tickv, textpos=0, minor=0)
c22=colorbar(target=im22, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c01.tickv, textpos=0, minor=0)
c23=colorbar(target=im23, orient=0, pos=[0, -0.07, 1, -0.02], /relative, $
  tickv=c03.tickv, textpos=0, minor=0)

t23 = text(mean(pos223[[0, 2]]), mean(pos222[[1, 3]]), /dev, $
    '$\sigma$ = '+string(stddev(diff, /nan), f='(f4.1)'), $
    align=0.5, vertical_align=0.5, $
    font_size=11, font_name='malgun gothic')


cd, '/hae/homedata/khcho/tomography-1906'
w01.save, 'fig2.pdf', resol=300, /bitmap, $
            page_size=w01.dimen/1d2, width=w01.dimen[0]/1d2
end

