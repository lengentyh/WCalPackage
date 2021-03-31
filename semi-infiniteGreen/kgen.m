clear all;

nkp=1;
kc=100;
k1=[-0.5  0.0  0.0   0.5  0.0  0.0];
% k2=[0.5  0.0  0.0   0.5  0.5  0.0];
% k3=[0.5  0.5  0.0   0.0  0.0  0.0];
% km=[k1;k2;k3];
km=k1;

kvector=[];
for ikk=1:nkp
    k=km(ikk,:);
    kx=linspace(k(1),k(4),kc);
    ky=linspace(k(2),k(5),kc);
    kz=linspace(k(3),k(6),kc);
    kvec=[kx' ky' kz'];
    kvector=[kvector; kvec];
end

kpath=kvector;
save kpath_slab.mat kpath
