odd_psi = fltarr(nx,ny,nz,ns)
even_psi = fltarr(nx,ny,nz,ns)
fft_odd = complexarr(nx,ny,ns)
fft_even = complexarr(nx,ny,ns)

for it=0,ns-1 do begin 
for ix=0,nx-1 do begin
odd_psi(ix,*,*,it)=(psi(nx-1-ix,*,*,it)-psi(ix,*,*,it))/2.
even_psi(ix,*,*,it)=(psi(nx-1-ix,*,*,it)+psi(ix,*,*,it))/2.
endfor

for ix=0,nx-1 do begin
fft_odd(ix,*,it) = fft(odd_psi(ix,*,0,it))
fft_even(ix,*,it) = fft(even_psi(ix,*,0,it))
endfor

endfor


end
