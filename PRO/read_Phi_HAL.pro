close,1

device,retain=2

nmfc='Phi.dat'

print,nmfc


time = 0.0d0
nx = 0L
ny = 0L
nz = 0L

stime=string(3)
openr,1,nmfc,/f77_unformatted
readu,1,time,nx,ny,nz
close,1

x = fltarr(nx)
work_gr = fltarr(1)

openr,3,'grid.dat'
for ix = 0, nx-1 do begin

readf,3,work_gr
x(ix) = work_gr(0)
endfor

close,3

print,'# time steps =  '
read,ns

print,'YL (-pi*yl<y<pi*yl)'
read,yl

print,'ZL (-pi*zl<y<pi*zl)'
read,zl

print,'# processors z =  '
read,nprocz
print,'# processors y =  '
read,nprocy

print,'time,nx,ny,nz,ns: '
print,time,nx,ny,nz,ns 

phi = fltarr(nx,ny,nz,ns)
vy = fltarr(nx,ny,nz,ns)
t = fltarr(ns)

y = (2.0 * !pi ) * yl * findgen(ny) / (ny) - yl * !pi 
z = (2.0 * !pi ) * zl * findgen(nz) / (nz) - zl * !pi

nzl = nz/nprocz
nyl = ny/nprocy

print,nzl

work = dblarr(nx,nyl,nzl)
my_t = fltarr(ns)

lx = 0L
ly = 0L
lz = 0L

openr,1,nmfc,/f77_unformatted

for it = 0,ns-1 do begin

readu,1,time,lx,ly,lz
print, time, lx, ly, lz

t(it) = time

for iprocy = 0,nprocy-1 do begin
for iprocz = 0,nprocz-1 do begin

readu,1,work

for iy = 0, ny/nprocy-1 do begin
for iz = 0, nz/nprocz-1 do begin
phi(*,iy + iprocy * nyl,iz + iprocz * nzl,it) = work(*,iy,iz)
vy(*,iy + iprocy * nyl,iz + iprocz * nzl,it) = deriv(x,phi(*,iy + iprocy * nyl,iz + iprocz * nzl,it))
endfor
endfor

endfor
endfor


endfor

close,1

end
