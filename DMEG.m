function [OBJ,W,pre_S1] = DMEG(X, Y, d,h,lambda ,maxiter)
% Initial projection matrix W
% X: data set, every colmun is a sample
% Y: label vector
% d: final dimension
% h: k nearest neighborhood
% Writen by Zheng Wang, email: zhengwangml@gmail.com
% S2 = cell(0);
[m,~] = size(X);
c = unique(Y);
n = length(Y);
H = eye(n) - 1/n*ones(n);
St = X*H*X'+ eps*eye(m);
invSt = inv(St);
pre_S = [];
Obj = 0;
OBJ = [];
% initual similarity matrix S
for i = 1 : length(c)
    Xc{i} = X(:,Y==i); 
    Xi = Xc{i};
    nc(i) = size(Xc{i},2);  
    distXi = L2_distance_1(Xi,Xi);
    [~, idx] = sort(distXi,2);
    S0{i} = construct_S0( idx, h, nc(i));   
    pre_S = blkdiag(pre_S,S0{i});      
end
pre_S1 = pre_S;
% Iterative calculate the projection W with S;   
OBJ = zeros(1,maxiter);

for iter = 1:maxiter
    S1 = [];
    % Update W    
    S  = (pre_S + pre_S')/2;
    D_s = diag(sum(S));
    L_w = D_s - S;
    Sw = X*L_w*X';
    P = invSt*Sw;  
    [W,~,~] = eig1(P, d, 0, 0);
    W = orth(W);
    WW = (W'*St*W)^(-0.5)*W';
    
    % Update S  
    for i=1:length(c)
    Xc{i} = X(:,Y==i); 
    nc(i) = size(Xc{i},2);
    Xi = Xc{i};
    ni = nc(i);
    distXi1 = L2_distance_1(WW*Xi,WW*Xi);
    [~, idx1] = sort(distXi1,2);
    SS{i} = construct_S6(distXi1, idx1, h, lambda, ni);
    obj(i) = sum(sum(distXi1.*SS{i})) + lambda * sum(SS{i}(SS{i}~=0).*log(SS{i}(SS{i}~=0)));
    S1 = blkdiag(S1,SS{i});
    end
    pre_S = S1;
    OBJ(iter) = sum(obj);
% %   Calculate objective function value
%     for i=1:length(c)
%     Xc{i} = X(:,Y==i); 
%     nc(i) = size(Xc{i},2);
%     Xi = Xc{i};
%     ni = nc(i);
%     distXi1 = L2_distance_1(WW1*Xi,WW1*Xi);
%     obj(i) = sum(sum(distXi1.*SS{i})) + lambda * sum(SS{i}(SS{i}~=0).*log(SS{i}(SS{i}~=0)));
% %   obj(i) = sum(sum(distXi.*SS{i}));
%     end
%     OBJ(iter) = sum(obj);
%     S2 = [S2,S1];
%     if 1 < iter && OBJ(iSter-1) < OBJ(iter) && iter < maxiter
%         disp(['Error: Algorithm does not reach convergence!!!','this is ', num2str(iter) ,' iteration'])
%     end
end 






