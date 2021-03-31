clear

%maxNumCompThreads(4);

load H_brick.mat
load T_brick.mat
load kpath_slab.mat

tic; % start clock

% ================================================================== %
%-------- input ---------%
Ni=10;                    % # of interation (layer); layer=2^Ni
EF=-3.7016;
w=EF+linspace(-0.5, 0.5,100);    % energy window
nn=0.01;                  % energy broadening
%EF=0;

%------------------------%

%kpath=kpath(:,1:2);      % path in k-space
nk=length(kpath);        % # of k
norb=ftn58_hop(1,1);     % # of orbital

%--- Choose the PDOS orbitals (OPTIONAL) --- %
PCD_orb        = 1:norb;
v              = ones(length(PCD_orb),1);
PCD_exp_sparse = sparse(PCD_orb,PCD_orb,v,norb,norb);
PCD_exp        = full(PCD_exp_sparse);

% =================================================================== %
% ------- Construct Matrix -------%
fprintf('Constructing matrix \n');
for ik=1:length(kpath)
% parfor ik=1:length(kpath)
    fprintf('matrix total k=%d, now=%d \n',length(kpath),ik);
    kcolumnvec=kpath(ik,:)'*2*pi;
    
    % ----- bulk part -----%
    % ----- transfor ftn58 to matirx (element) ----%
    
    % ---- diagonal element H ----%
    Hk{ik}=full(sparse(ftn58_brick(2:end,2),...
        ftn58_brick(2:end,3),...
        exp(i*ftn58_brick(2:end,5:7)*kcolumnvec).*ftn58_brick(2:end,4),...
        norb,norb));
    Hk{ik}=(Hk{ik}+Hk{ik}')/2;   % keep Hermination
    
    %     % --- off-diagonal term T ---%
    Tk{ik}=full(sparse(ftn58_hop(2:end,2),...
        ftn58_hop(2:end,3),...
        exp(i*ftn58_hop(2:end,5:7)*kcolumnvec).*ftn58_hop(2:end,4),...
        norb,norb));
    
    %     % ---- surface part ----%
    %     % ---- diagonal element H ----%
    Hs0{ik}=full(sparse(ftn58_brick(2:end,2),...
        ftn58_brick(2:end,3),...
        exp(i*ftn58_brick(2:end,5:7)*kcolumnvec).*ftn58_brick(2:end,4),...
        norb,norb));
    Hs0{ik}=(Hs0{ik}+Hs0{ik}')/2;
    % --- off-diagonal term T ---%
    Ts0{ik}=full(sparse(ftn58_hop(2:end,2),...
        ftn58_hop(2:end,3),...
        exp(i*ftn58_hop(2:end,5:7)*kcolumnvec).*ftn58_hop(2:end,4),...
        norb,norb));
    
end

%return

% Sx
dnup=eye(norb/2);
updn=dnup;
sx=[zeros(norb/2), dnup; updn, zeros(norb/2)];
% Sy
dnup=-i*eye(norb/2);
updn=-dnup;
sy=[zeros(norb/2), dnup; updn, zeros(norb/2)];
% Sz
upup=eye(norb/2);
dndn=-upup;
sz=[upup, zeros(norb/2); zeros(norb/2), dndn];


% ---Memory space ----%
II=eye(norb);

Ab=zeros(length(w),length(kpath));
As=zeros(length(w),length(kpath));
%As_s=zeros(length(w),length(kpath));
As_sx=zeros(length(w),length(kpath));
As_sz=zeros(length(w),length(kpath));
As_sy=zeros(length(w),length(kpath));
%---------------------%

ef=EF;

%------- Green function and spectral weight ----%
for iw=1:length(w)    % scan energy
% parfor iw=1:length(w)    % scan energy
    fprintf('total w=%d, now=%d \n',length(w),iw);
    runw=w(iw);
    Abt=Ab(iw,:);
    Ast=As(iw,:);
    %    Ast_s=As_s(iw,:);
    Ast_sx=As_sx(iw,:);
    Ast_sy=As_sy(iw,:);
    Ast_sz=As_sz(iw,:);
    
    for ik=1:length(kpath)    % scan k
        E=(runw+1i*nn)*II;    % w in geen function
        H=full(Hk{ik});
        T=full(Tk{ik});
        H0=full(Hs0{ik});
        T0=full(Ts0{ik});
        [Hs,Hb]=iterationH(E,H,T,Ni);  % PRB 31, 5166 eq.10
        
        Gb=inv(E-Hb);  % bulk  % PRB 31, 5166 eq.5
        Gs=inv(E-Hs);  % surface
        
        Gs0=inv(E-H0-T0*Gs*T0');    % cover surface % PRB 31, 5166 eq.13
        
        Abt(ik)=-imag(trace(Gb*PCD_exp));   % spectral weight of bulk
        Ast(ik)=-imag(trace(Gs0*PCD_exp));  % spectral weight of surface (sun over all of total orbital)
        %        Ast(ik)=-imag(trace(Gs0(49:52,49:52))+trace(Gs0(121:124,121:124)));
        %                        trace(Gs0(129:144,129:144))+trace(Gs0(321:336,321:336)));
        
        Gs0sx=sx*Gs0;
        Gs0sy=sy*Gs0;
        Gs0sz=sz*Gs0;
        Ast_sx(ik)=-imag(trace(Gs0sx*PCD_exp));
        Ast_sy(ik)=-imag(trace(Gs0sy*PCD_exp));
        Ast_sz(ik)=-imag(trace(Gs0sz*PCD_exp));
        
        %        Ast_sx(ik)=-imag(trace(Gs0sx(49:52,49:52)+Gs0sx(121:124,121:124)));
        %        Ast_sy(ik)=-imag(trace(Gs0sy(49:52,49:52)+Gs0sy(121:124,121:124)));
        %        Ast_sz(ik)=-imag(trace(Gs0sz(49:52,49:52)+Gs0sz(121:124,121:124)));
        
    end
    %Ab(iw,:)=Abt;
    
    Ab(iw,:)=Abt;
    As(iw,:)=Ast;
    %    As_s(iw,:)=Ast_s;
    As_sx(iw,:)=Ast_sx;
    As_sy(iw,:)=Ast_sy;
    As_sz(iw,:)=Ast_sz;
    
end

%At=As+Ab;

%w=w-ef;

% ----------------------------------- %

save ek_density.mat  w nk ef Ab As
save ek_spin.mat  As_sx As_sy As_sz
% save G_matrix.mat -v7.3 Gb

toc; % stop clock
Computation_time=num2str(toc); % computation time


return
%%
load ek_density.mat
wef=w-ef;
AA=As; % Ab:bulk As:surface
figure; hold on
kd=linspace(0,nk,nk);
nk=length(kd);
xlabel('K_x')
wef=w-ef;
ymin=-0.5; ymax=0.5;
set(gca,'XTick',[],'XTickLabel',{''})


pcolor(kd,wef,AA/pi),shading interp

axis square,box on ,set(gca,'linewidth',2,'FontSize',16)
plot([kd(1) kd(end)],[ 0 0],'k-.')

axis([0 kd(end) ymin ymax])
ylabel('Energy (eV)')
