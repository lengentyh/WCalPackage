function A=bandplot(ftn58, kpt)

%global norb Hfix icross ftn58 subindex

%%% hand make input %%%

%kpt = buildkpts('YBCO123','subdir10');

%kpt = buildkpts('YBCO123_supercell_2','subdir');
%kpt = buildkpts('YBCO123_afm_5','subdir');

kpoints = kpt(:,3:5);
norb=ftn58(1,1);
kmn=size(kpoints,1);
A=zeros(kmn*norb,2);
for ik=1:kmn
    %A(ik*norb-norb+1:ik*norb,1)=ik;
    A(ik*norb-norb+1:ik*norb,1)=kpt(ik,2);
%     kx=klist(ik,1)/klist(ik,4)*2*pi;
%     ky=klist(ik,2)/klist(ik,4)*2*pi;
%     kz=klist(ik,3)/klist(ik,4)*2*pi;
    A(ik*norb-norb+1:ik*norb,2)=eigband(ftn58,kpoints(ik,:)*2*pi);
end
return
