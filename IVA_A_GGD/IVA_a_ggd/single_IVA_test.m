
% For a general description of the IVA-A-GGD algorithm and its relationship with others,
% see http://mlsp.umbc.edu/jointBSS_introduction.html

clear all
close all
clc

N=9;
K=3;
T=4096;

S=zeros(N,T,K);

for n=1:N
    beta(n) = rand*(4-0.25)+0.25;
    rho(n) = rand*(0.7-0.5)+0.5; %Correlation within the SCV
    Sigma = cov_sigma(K, rho(n)); % scatter matirx
    %U = rand(K);
    %Sigma = U*U';
    S(n,:,:)=generate_MGGD(T, K, Sigma, beta(n))';
    %S(n,:,:)=randmv_laplace(K,T)';
    %S(n,:,:) = (Z(1:K/2,:) + 1i*Z(K/2+1:end,:)).';
end
A=randn(N,N,K);% + 1j*randn(N,N,K/2);

for k=1:K
    %A(:,:,k)=vecnorm(A(:,:,k))';
    X(:,:,k)=A(:,:,k)*S(:,:,k);
end

%%
for kk=1:K
    wnit(:,:,kk) = rand(N,N);
end

tic
[W1, cost1, shapMoM, isi1, iterMoM] = iva_a_ggd_decp_RA_FP(X,'A', A,'initW',wnit);
toc

tic
[W3,cost3,shapeParam3,isi3] = iva_a_ggd_decp_RA_FP_P(X,'A',A,'initW',wnit);
toc
