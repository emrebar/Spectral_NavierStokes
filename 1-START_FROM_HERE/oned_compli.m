 %% BD 3
        clc;clear all;close all;
        Ny=256;
        kx = 4;
        [D,y] = cheb(Ny);           %cos(pi*(0:Ny)/Ny);
        yi = linspace(-1,1,1000);
        D2 = D*D;
        ufun = @(y,t) sin(y.^3+t)+5*t;
        ufun_y = @(y,t) 3*y.*cos(y.^3+t);
        ufun_t = @(y,t) cos(y.^3+t)+5;
        ufun_yy = @(y,t) 3.*cos(y.^3 + t) - 9*y.^3.*sin(y.^3 + t);
        bfun = @(y,t) ufun_t(y,t) - ufun_yy(y,t) + kx^2*ufun(y,t) ;
        u = ufun(y,0);
        I=eye(Ny+1);
        p=1e1;     % store
        step=5e5;
        dt= 1e-6;
        t=(0:step)*dt;
        U=zeros(step,Ny+1);
        uold=ufun(y,-dt);
        uoldd=ufun(y,-2*dt);
        
        A=(11/(6*dt))*I-D2+kx^2*I;
        
        for j=1:step;
            c=(3/dt)*u-(3/(2*dt))*uold+1/(3*dt)*uoldd+bfun(y,t(j)+dt);
            unew=A\c;
            unew(1)=ufun(y(1),t(j));    unew(end)=ufun(y(end),t(j));
            uoldd=uold;uold=u;u=unew;
            
            %if mod(step,p)==0;     U(j,:)=u; end
            U(j,:)=u; 
        end
%         bd31dlong=U;
%         save bd31dlong bd31dlong
% end


%% Visu.
% load rk4 bd3 bd2
figure(1),clf(1)
for ii=1:100:length(t)-1
    plot(U(ii,:),y,'b.-',ufun(yi,t(ii)),yi,'r-')
    title(sprintf('time=%-4.3e',t(ii))),drawnow,shg
end