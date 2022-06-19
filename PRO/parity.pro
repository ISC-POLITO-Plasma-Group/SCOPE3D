odd_phi = fltarr(nx,ny,nz,ns)
even_phi = fltarr(nx,ny,nz,ns)
fft_odd = complexarr(nx,ny,ns)
fft_even = complexarr(nx,ny,ns)

for it=0,ns-1 do begin 
for ix=0,nx-1 do begin
odd_phi(ix,*,*,it)=(phi(nx-1-ix,*,*,it)-phi(ix,*,*,it))/2.
even_phi(ix,*,*,it)=(phi(nx-1-ix,*,*,it)+phi(ix,*,*,it))/2.
endfor

for ix=0,nx-1 do begin
fft_odd(ix,*,it) = fft(odd_phi(ix,*,0,it))
fft_even(ix,*,it) = fft(even_phi(ix,*,0,it))
endfor

endfor


end
