function ph=aapx(uh,vh) %anti-aliased product, 
% in: uh,vh from fft with n samples
[ny nx]=size(uh);m=nx*3/2;
uhp=[uh(:,1:nx/2) zeros(ny,(m-nx)) uh(:,nx/2+1:nx)]; % pad uhat with zeros 
vhp=[vh(:,1:nx/2) zeros(ny,(m-nx)) vh(:,nx/2+1:nx)]; % pad vhat with zeros 
up=ifft(uhp,[],2); vp=ifft(vhp,[],2); w=up.*vp; wh=fft(w,[],2);
ph=1.5*[wh(:,1:nx/2) wh(:,m-nx/2+1:m)]; % extract F-coefficients 
end 