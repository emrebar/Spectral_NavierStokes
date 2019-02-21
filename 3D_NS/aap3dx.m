function ph=aap3dx(uh,vh) %anti-aliased product, 
% in: uh,vh from fft with n samples
[ny nx nz]=size(uh);m=nx*3/2;mm=nz*3/2;
uhp=[uh(:,1:nx/2,:) zeros(ny,(m-nx),nz) uh(:,nx/2+1:nx,:)]; % pad uhat with zeros 
vhp=[vh(:,1:nx/2,:) zeros(ny,(m-nx),nz) vh(:,nx/2+1:nx,:)]; % pad vhat with zeros 
uhhp=shiftdim(uhp,1);vhhp=shiftdim(vhp,1);
uhpn=[uhhp(:,1:nz/2,:) zeros(m,(mm-nz),ny) uhhp(:,nz/2+1:nz,:)]; % pad uhat with zeros 
vhpn=[vhhp(:,1:nz/2,:) zeros(m,(mm-nz),ny) vhhp(:,nz/2+1:nz,:)]; % pad uhat with zeros 
uhpn=shiftdim(uhpn,2);vhpn=shiftdim(vhpn,2);
up=ifft(ifft(uhpn,[],3),[],2); vp=ifft(ifft(vhpn,[],3),[],2); conv=up.*vp; convh=fft(fft(conv,[],3),[],2);
ph=1.5*[convh(:,1:nx/2,:) convh(:,m-nx/2+1:m,:)]; % extract F-coefficients 
ph=shiftdim(ph,1);
ph=1.5*[ph(:,1:nz/2,:) ph(:,mm-nz/2+1:mm,:)]; % extract F-coefficients 
ph=shiftdim(ph,2);
end 