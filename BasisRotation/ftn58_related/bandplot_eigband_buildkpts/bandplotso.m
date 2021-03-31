function A=bandplotso(ftn58, kpoints, ftn58SO)

%global norb Hfix icross ftn58 subindex


norb=ftn58(1);
kmn=size(kpoints,1);
A=zeros(kmn*norb,2);
for ik=1:kmn
    A(ik*norb-norb+1:ik*norb,1)=ik;
%     kx=klist(ik,1)/klist(ik,4)*2*pi;
%     ky=klist(ik,2)/klist(ik,4)*2*pi;
%     kz=klist(ik,3)/klist(ik,4)*2*pi;
    A(ik*norb-norb+1:ik*norb,2)=eigbandso(ftn58,kpoints(ik,:)*2*pi),ftn58SO);% For (x,y,z) system
    %A(ik*norb-norb+1:ik*norb,2)=eigband(ftn58,kpoints(ik,:)*2*pi); % For (a1,a2,a3) system
end
return
