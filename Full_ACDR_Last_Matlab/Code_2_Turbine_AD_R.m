load preprocess;load tjaturbine.mat
%% Turbine airfoil & Radius
SampleTurbineR=30.46;
%% Smear Global Parameters
epi=1.1;ep=epi*dx;
%% Find the CENTERS within rotor disk
coord=[X(:) Y(:) Z(:)];
xcenter=find(coord(:,1)==x(ceil((Nx/5))));
[YY]=coord(xcenter,2);[ZZ]=coord(xcenter,3);
YY=YY-ygl(Ny/2);ZZ=ZZ-z(Nz/2);
hhh=zeros(length(YY));
for i=1:size(YY(:));zz=ZZ(i);yy=YY(i);
    rad=norm([zz yy]);if rad<0.5;hhh(i)=1;end;
end
circle=find(hhh);
CENTERS=[x(ceil(Nx/5))*ones(size(circle)) YY(circle)+ygl(Ny/2) ZZ(circle)+z(Nz/2)];
%% Find the indices where the turbine is located and store in turbindex
turbindex=zeros(length(CENTERS),3);
for ll=1:length(CENTERS); 
     y0=CENTERS(ll,2);z0=CENTERS(ll,3);
    yncent(ll)=y0-ygl((Ny/2));zncent(ll)=z0-z((Nz/2));
    radinside(ll)=norm([yncent(ll),zncent(ll)]);
turbindex(ll,:)=[find(x==CENTERS(ll,1))   find(ygl==CENTERS(ll,2))   find(z==CENTERS(ll,3))];
end;

%% Pick correct Airfoil Profiles for correct thickness store in profcatalog
for kkk=1:length(CENTERS);
   localthick(kkk)=interp1(tempbladerad/SampleTurbineR,[toverc(1) toverc],(radinside(kkk)/0.5));
end
simulationturbine=[12 15 18 21 24]./100;
for kkk=1:length(CENTERS);
aaa=abs(simulationturbine-localthick(kkk));
[profile]=find(aaa==min(aaa));
profcatalog(kkk,:)=[kkk simulationturbine(profile)];
end



save ('turbines','CENTERS','turbindex','profcatalog','ep','radinside','yncent','zncent','SampleTurbineR')