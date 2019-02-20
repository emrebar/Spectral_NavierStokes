function [turbxx,turbyy,turbzz]=fxftheta(u,v,w)

global CENTERS x ygl z X Y Z Nx Ny j Nz turbindex omega dx ep f_eps f1d vol
global alfclcd12 alfclcd15 alfclcd18 alfclcd21 alfclcd24 bladerad chordtjae elementarea tempbladerad thicktja twisttja toverc
global aind apind radinside SampleTurbineR profcatalog yncent zncent f_eps1 f_eps2

turbxx=0;turbyy=0;turbzz=0;

for kkk=1:length(CENTERS)
    
    % assign computed velocities on the domain to Vx(kkk),Vy,Vz
    Vx(kkk)=-u(turbindex(kkk,2),turbindex(kkk,1),turbindex(kkk,3));
    Vy=-v(turbindex(kkk,2),turbindex(kkk,1),turbindex(kkk,3));
    Vz=-w(turbindex(kkk,2),turbindex(kkk,1),turbindex(kkk,3));
    % find the angle to project Vtheta(kkk)
    beta=atand(Vy/Vz);
    % project to Vtheta(kkk)
    Vthetaprojected=dot([0 -sind(beta) cosd(beta)],[0 Vy Vz]);
    
    %     if isnan(Vtheta(kkk)projected)==1; Vtheta(kkk)projected=1.5;end
    
    % take into accoun the rotation of the turbine and find the real Vtheta(kkk)
    Vtheta(kkk)=omega*(radinside(kkk)/0.5)*SampleTurbineR+Vthetaprojected;
    if j<300;    Vx(kkk)=1-aind;   Vtheta(kkk)=(omega*(radinside(kkk))/0.5*30.46)*(1+apind);end
    
    % find relative velocity for lift & drag calculation
    Vrel=norm([Vx(kkk),Vtheta(kkk)]);
    
    
    
    % find inflow angle
    phi=atand(Vx(kkk)/Vtheta(kkk));
    % find localpitch(kkk)
    localpitch(kkk)=interp1(bladerad/SampleTurbineR,twisttja,(radinside(kkk)/0.5));
    % angle of attack= inflow angle - (localpitch(kkk)+twist)
    aoa(kkk)=phi-localpitch(kkk);
    %     if aoa(kkk)>10;aoa(kkk)=10;end
    % find which thickness airfoil is used in that section
    if profcatalog(kkk,2)>1;profcatalog(kkk,2)=0.12;end
    localairfoilprof=eval(['alfclcd',num2str(profcatalog(kkk,2)*100)]);
    % find Cl & Cd
    Cllocal(kkk)=interp1(localairfoilprof(:,1),localairfoilprof(:,2),aoa(kkk));
    Cdlocal(kkk)=interp1(localairfoilprof(:,1),localairfoilprof(:,3),aoa(kkk));
    % find local-area
    localarea(kkk)=interp1(bladerad/SampleTurbineR,elementarea,(radinside(kkk)/0.5));
    
    % calculate lift & drag
    lift=((Vrel^2)/2)*localarea(kkk)*Cllocal(kkk);
    drag=((Vrel^2)/2)*localarea(kkk)*Cdlocal(kkk);
    
    % project the force on the x direction
    Fx(kkk)=(lift*cosd(phi)+drag*sind(phi))/(2*pi*(radinside(kkk)/0.5)*SampleTurbineR);
    Ftheta(kkk)=(lift*sind(phi)-drag*cosd(phi))/(2*pi*(radinside(kkk)/0.5)*SampleTurbineR);
    
    %  add nacelle- adjust the area and Cd accordingly.
    if radinside(kkk)<0.10;Fx(kkk)=0.1*0.1*0.1/2/vol*dx;
        
    end
    
    %% Smeared forces
    if kkk<75;
    turbxx=f_eps1{kkk}*Fx(kkk)/dx+turbxx;
    if isnan(Ftheta(kkk))==1;Ftheta(kkk)=Ftheta(kkk-1);end
    if radinside(kkk)==0;sortedrad=sort(radinside);radinside(kkk)=sortedrad(2); end

    turbyy=f_eps1{kkk}*(Ftheta  (kkk)/dx)*  (-zncent(kkk)/radinside(kkk))+turbyy;
    turbzz=f_eps1{kkk}*(Ftheta  (kkk)/dx)*  (yncent(kkk)/radinside(kkk))+turbzz;
    
    else
        turbxx=f_eps2{kkk-74}*Fx(kkk)/dx+turbxx;
    if isnan(Ftheta(kkk))==1;Ftheta(kkk)=Ftheta(kkk-1);end
    if radinside(kkk)==0;sortedrad=sort(radinside);radinside(kkk)=sortedrad(2); end

    turbyy=f_eps2{kkk-74}*(Ftheta  (kkk)/dx)*  (-zncent(kkk)/radinside(kkk))+turbyy;
    turbzz=f_eps2{kkk-74}*(Ftheta  (kkk)/dx)*  (yncent(kkk)/radinside(kkk))+turbzz; 
        
    end
end
%%
% turbxx=0;for kkk=1:length(CENTERS);turbxx=turbx{kkk}+turbxx;end
% turbyy=0;for kkk=1:length(CENTERS);turbyy=turby{kkk}+turbyy;end
% turbzz=0;for kkk=1:length(CENTERS);turbzz=turbz{kkk}+turbzz;end
% [sum(abs(turbyy(:))) sum(abs(turbzz(:))) sum(abs(turbxx(:)))];


 end