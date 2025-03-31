cd,'/data/home/chokh/IDL_default/idl_lib/tomography'
file = FILEPATH('head.dat', SUBDIRECTORY = ['examples', 'data'])
data = READ_BINARY(file, DATA_DIMS=[80, 100, 57])

sz = 128.
head3d = congrid(data, sz, sz, sz)

xp = findgen(sz)
xc = (sz-1.)*0.5
xxp = rebin(xp, sz, sz, sz)
yyp = rebin(reform(xp, 1, sz), sz, sz, sz)
outside = where((xxp-xc)^2.+(yyp-xc)^2. gt (sz*0.5)^2.)
head3d[outside] = 0.

save, sz, head3d, filename='head3d.sav'

end