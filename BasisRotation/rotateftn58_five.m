function [new_ftn58sparse,trot]=rotateftn58_five(ftn58sparse,a,b,c,orbitype,orbital)
    % ftn58 should be in the Workspace 
    % a,b,c are Euler angles with rotation Rz(a)*Ry(b)*Rz(c)
    % orbitype = {'p','d','f'}
    % orbital = [1,2,3], i.e., [1st, 2nd, 3rd] orbital on Wannier90
    % These orbitals should be in a continuous order on Wannier90
    % And these orbitals should be ordered in the following sequence 
    % p:pz,px,py
    % d:dz2, dxz, dyz, dx2-y2, dxy
    % f:fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx(x2-3y2), fy(3x2-y2)
    

   % initialization 
   new_ftn58sparse = ftn58sparse;
   ij 		   = ftn58sparse.ij;
   tt		   = ftn58sparse.tt;
   dd		   = ftn58sparse.dd;
   norb            = ftn58sparse.norb;
   nbond           = size(ij,1);
    
    % add column in ftn58 to classify translation
    trsl=unique(dd,'rows');
    [logn,ordn] = ismember(dd,trsl,'rows');
    
    % pull out and construct hopping matrix
    nlen = length(trsl);
    t = zeros(norb,norb,nlen);
    for j=1:nbond  
        m=ij(j,1);
        k=ij(j,2);
        n=ordn(j);
        t(m,k,n)=tt(j);
    end
    
    % rotate (t) for each cell
    D = genWD(a,b,c,orbitype);
    I1 = eye(min(orbital)-1);
    I2 = eye(norb-max(orbital));
    ROT = blkdiag(I1,D,I2);
    trot = zeros(norb,norb,nlen);
    for j=1:nlen
        trot(:,:,j)=ctranspose(ROT)*t(:,:,j)*ROT;
    end
    
    % put back trot into ftn58
    for j=1:nbond
        new_ftn58sparse.tt(j,1) = trot(ij(j,1),ij(j,2),ordn(j));
    end

end
