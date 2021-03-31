clear all;

load ftn58soc_full.mat
ftn58=ftn58soc_full;
%% <!> ftn58 must be hermitian !!! 
%% Since the following codes doesn't make T_brick to be hermitian !!! 

% ----- Input parameter -----%
surface=2; % surface x=1, y=2, z=3
hop_d=1;   % hopping direction. positive=1(coord small surface), negeative=2;
%e.g hop_d = 1 means that the lower atoms hop to highet atoms , s.t the
%R>0, and using these T to calculate semi-infinite Green , we have the 
%semi-infinite Green function which terminated at the lower surface
%------------- if hop_d=2(semi-infinite at this surface)
%   BulkBulkBulk
%   BulkBulkBulk
%   BulkBulkBulk
%   BulkBulkBulk
%------------- if hop_d=1(semi-infinite at this surface)
%
%
%----------------------------%

if surface==1;
    surf_d=5;  % ftn58 x
end
if surface==2;
    surf_d=6;  % ftn58 y
end
if surface==3;
    surf_d=7;  % ftn58 z
end

% longest hopping in bulk
%long_hop_p=max(max(ftn58(:,surf_d))); % positive hop 
%long_hop_n=min(min(ftn58(:,surf_d))); % negative hop

%%% --- to obtain the propagator (its ftn58)-----------------------------
%%% --- for calculation of the renormalization --------------------------
%%% --- the longest hopping !!! -----------------------------------------
long_hop_p=1; % positive hop 
long_hop_n=-1; % negative hop

% brick 
%brick_ind=find(ftn58(:,surf_d)<long_hop_p & ftn58(:,surf_d)>long_hop_n); % search element of H1 in PRB 31, 5166 eq.12

% the ftn58 of the 
if hop_d==1;
    brick_ind = find(ftn58(:,surf_d)<long_hop_p & ftn58(:,surf_d)>=0);
else
    brick_ind = find(ftn58(:,surf_d)>long_hop_n & ftn58(:,surf_d)<=0);
end

ftn58_brick=ftn58(brick_ind(:),:);
%n_element_brick=size(brick_ind);         % # of elements
%ftn58_brick=zeros(n_element_brick(1),7); % memory
%for iib=1:n_element_brick(1)             % element of brick (ftn58 fromat)
%    ftn58_brick(iib,:)=ftn58(brick_ind(iib),:);
%end

% positive hopping (brick to brick)
if hop_d==1;
    hop_ind_p=find(ftn58(:,surf_d)>0);    % search element of T1 in PRB 31, 5166 eq.12
    ftn58_hop_p=ftn58(hop_ind_p(:),:);
    %n_element_hop_p=size(hop_ind_p);      % # of elements 
    %ftn58_hop_p=zeros(n_element_hop_p(1),7); % memory 
    %for iih=1:n_element_hop_p(1)             % element of hopping (ftn58 fromat)  
    %    ftn58_hop_p(iih,:)=ftn58(hop_ind_p(iih),:);
    %end
    ftn58_hop=[ftn58(1,:); ftn58_hop_p];
else
% negative hopping (brick to brick)
    hop_ind_n=find(ftn58(:,surf_d)<0);    
    ftn58_hop_n=ftn58(hop_ind_n(:),:);
    %n_element_hop_n=size(hop_ind_n);
    %ftn58_hop_n=zeros(n_element_hop_n(1),7);
    %for iih=1:n_element_hop_n(1)
    %    ftn58_hop_n(iih,:)=ftn58(hop_ind_n(iih),:);
    %end
    ftn58_hop=[ftn58(1,:); ftn58_hop_n];
end

save H_brick.mat ftn58_brick
save T_brick.mat ftn58_hop




