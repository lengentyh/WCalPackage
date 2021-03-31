function [Lz,Lp]=genL(orbit)
    % generate Lz, L+ matrix for p, d, f orbit
    % orbit={'p','d','f'}; 
	% p:pz,px,py
	% d:dz2, dxz, dyz, dx2-y2, dxy
	% f:fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx(x2-3y2), fy(3x2-y2)
    Lz=[];
    Lp=[];
    if any(strcmp(orbit,'p'))
        Lzt=zeros(3);
        Lzt(2,3)=-1i;
        Lzt=Lzt+Lzt';
        Lz=cat(1,Lz,{Lzt});
        Lpt=zeros(3);
        Lpt(1,2:3)=[1 1i];
        Lpt(2:3,1)=[-1 -1i];
        Lp=cat(1,Lp,{Lpt});
    end
    if any(strcmp(orbit,'d'))
        Lzt=zeros(5);
        Lzt(2,3)=-1i;
        Lzt(4,5)=-2i;
        Lzt=Lzt+Lzt';
        Lz=cat(1,Lz,{Lzt});
        Lpt=zeros(5);
        Lpt(1,2:3)=sqrt(3)*[1 1i];
        Lpt(2:3,1)=sqrt(3)*[-1 -1i];
        Lpt(2:3,4:5)=[1 1i;-1i 1];
        Lpt(4:5,2:3)=[-1 1i;-1i -1];
        Lp=cat(1,Lp,{Lpt});
    end
	if any(strcmp(orbit,'f'))
		Lzt=[0 -1i 0 -2i 0 -3i];
		Lzt=diag(Lzt,1);
		Lzt=Lzt+Lzt';
		Lz=cat(1,Lz,{Lzt});
        Lpt=zeros(7);
        Lpt(1,2:3)=[1 1i]*sqrt(6);
        Lpt(2:3,4:5)=[1 1i;-1i 1]*sqrt(5/2);
        Lpt(4:5,6:7)=[1 1i;-1i 1]*sqrt(3/2);
        Lpt=Lpt-Lpt.';
        Lp=cat(1,Lp,{Lpt});
	end
end