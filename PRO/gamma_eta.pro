close,1

device,retain=2

nmfc='Phi.dat'

print,nmfc


time = 0.0d0
nx = 0L
ny = 0L
nz = 0L
m = 12.0
toll = 0.00001

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
odd_phi = fltarr(nx,ny,nz,ns)
even_phi = fltarr(nx,ny,nz,ns)
fft_odd= complexarr(nx,ny,ns)
fft_even = complexarr(nx,ny,ns)
t = fltarr(ns)
a = fltarr(m,ns)
b= fltarr(m)
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

for ix=0,nx-1 do begin
odd_phi(ix,*,*,it)=(phi(nx-1-ix,*,*,it)-phi(ix,*,*,it))/2.
even_phi(ix,*,*,it)=(phi(nx-1-ix,*,*,it)+phi(ix,*,*,it))/2.
endfor

for ix=0,nx-1 do begin
fft_odd(ix,*,it)=fft(odd_phi(ix,*,0,it)) 
fft_even(ix,*,it)=fft(even_phi(ix,*,0,it))
endfor



endfor

for im=0,m-1 do begin
a(im,*) = deriv(t,alog(abs(fft_odd(nx/2,im,*))))
b(im) = a(im,ns-2)
endfor

plot,b(1:m-1)

print, 'growth rate'

print,b(1:m-1)

corr = fltarr(m)
difference = fltarr(m)
corr = [0,0.0150459,0.0225031,0.0240033,0.021275,0.0174198,0.0149267,0.0136614,0.0123065,0.0108228,0.00878358,0.00589085]
print, 'correct growth rate='
print, corr(1:m-1)
for im=1,m-1 do begin
difference(im) =b(im)-corr(im)

if (difference(im) gt toll) then begin
print, 'TEST FAILED'
endif else begin
print, 'TEST PASSED'
endelse

endfor
print,'difference='
print, difference(1:m-1)


close,1

end
