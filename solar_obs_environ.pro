path='/data/home/chokh/IDL_default/idl_lib/tomography'
restore, 'sinogram1.sav', /verbose
;stop
sz=size(sinogram)
xc=(sz[1]-1)*0.5
yc=(sz[1]-1)*0.5
xp=findgen(sz[1])-xc
yp=findgen(sz[1])-yc
xxp=rebin(xp, sz[1], sz[1])
yyp=rebin(transpose(yp), sz[1], sz[1])
data[where(sqrt(xxp^2.+yyp^2.) lt r_sun)] = 255.

w01 = window(dim=[8d2, 8d2])
im01 = image(data, pos=[0.2, 0.2, 0.8, 0.8], /current)
p01 = plot((xc-r_sun)*[1, 1], [0, sz[1]], over=im01, $
            '--3', color='orange')
p02 = plot((xc+r_sun)*[1, 1], [0, sz[1]], over=im01, $
            '--3', color='orange')            

a01 = arrow(xc*[1, 1], sz[1]+[5, 35], clip=0, /data, /current, $
            thick=2) 
t01 = text(xc, sz[1]+40, 'Observer', /data, clip=0, $
           font_size=15, font_style=1, font_name='malgun gothic', $
           align=0.5)

t02 = text(0.5*(xc-r_sun), 220, '(a)', /data, color='yellow',$
           font_size=15, font_style=1, font_name='malgun gothic', $
           align=0.5)

t03 = text(xc, 220, '(b)', /data, color='yellow',$
           font_size=15, font_style=1, font_name='malgun gothic', $
           align=0.5)
           
t04 = text(xc+0.5*(xc+r_sun), 220, '(c)', /data, color='yellow',$
           font_size=15, font_style=1, font_name='malgun gothic', $
           align=0.5)

t05 = text(xc, yc, 'Solar !cInterior', /data, $
           font_size=15, font_style=1, font_name='malgun gothic', $
           align=0.5, vertical_align=0.5)
           
t05 = text(xc, 30, '(d)', /data, color='yellow',$
           font_size=15, font_style=1, font_name='malgun gothic', $
           align=0.5)

w01.save, 'illustration.pdf', page_size=8.5*[1, 1], width=8.5
end