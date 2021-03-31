function D=genWD(a,b,c,orbit)
    % generate Wigner D matrice according to given Euler angles
    % rotations: Rz(a)*Ry(b)*Rz(c)
    % orbit={'p','d','f'};
    % p:pz,px,py
	% d:dz2, dxz, dyz, dx2-y2, dxy
	% f:fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx(x2-3y2), fy(3x2-y2)
   [Lz,Lp]=genL(orbit);
   Lz = cell2mat(Lz);
   Lp = cell2mat(Lp);
   Ln=ctranspose(Lp);
   D=expm(-1i*a*Lz)*expm(-1i*b*(1i*(Ln-Lp)/2))*expm(-1i*c*Lz);
end