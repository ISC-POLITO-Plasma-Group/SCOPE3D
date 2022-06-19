close,1

device,retain=2

nmfc='jre.dat'

print,nmfc


time = 0.0d0
nx = 0L
ny = 0L
nz = 0L

stime=string(3)
;openr,1,nmfc,/f77_unformatted,/Swap_If_little_Endian
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

print,'eq_L (psi_eq = -A*alog(cosh(x/eq_L)))'
read,eq_l

print,'A (-A*alog(cosh((x/eq_L)))'
read,eq_ampl

print,'asym (n/m)'
read,asym

print,'# processors z =  '
read,nprocz
print,'# processors y =  '
read,nprocy

print,'time,nx,ny,nz,ns: '
print,time,nx,ny,nz,ns 

jre = fltarr(nx,ny,nz,ns)
jre_tot = fltarr(nx,ny,nz,ns)
by = fltarr(nx,ny,nz,ns)
bx = fltarr(nx,ny,nz,ns)
;d2x = fltarr(nx,ny,nz,ns)
;d2y = fltarr(nx,ny,nz,ns)
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

;openr,1,nmfc,/f77_unformatted,/Swap_If_little_Endian
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
jre(*,iy + iprocy * nyl,iz + iprocz * nzl,it) = work(*,iy,iz)
jre_tot(*,iy + iprocy * nyl,iz + iprocz * nzl,it) = work(*,iy,iz)+eq_ampl/(cosh(x(*)/eq_l))+asym*yl/zl*x(*)
by(*,iy + iprocy * nyl,iz + iprocz * nzl,it) = deriv(x,jre_tot(*,iy + iprocy * nyl,iz + iprocz * nzl,it))
endfor
endfor

endfor
endfor

;for iz = 0, nz-1 do begin
;for iy = 0, ny-1 do begin
;d2x(*,iy,iz,it) = deriv(x,deriv(x,psi(*,iy,iz,it)))
;endfor
;for ix = 0, nx-1 do begin
;d2y(ix,*,iz,it) = deriv(y,deriv(y,psi(ix,*,iz,it)))
;endfor
;endfor

;for iz = 0, nz-1 do begin
;for ix = 0, nx-1 do begin
;bx(ix,*,iz,it) = deriv(y,psi(ix,*,iz,it))
;endfor
;endfor


endfor

close,1

end
