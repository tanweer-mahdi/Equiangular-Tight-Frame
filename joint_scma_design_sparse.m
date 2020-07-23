clc;
close all;
clear all;
% rng('shuffle');
%% This code develops Joint SCMA Codebook where the codewords has low mutual coherence %%

%% Parameter initialization
N = 24; %Number of resources assigned per user
M = 60; %codebook size
K = 2*N; %number of resources
overloading_factor = 1.5;
J = overloading_factor*K; %number of user
%Es = (N/2); %average symbol energy
Es = 1;
%% Mapping Matrix Initialization

V=zeros(K,N,J); % Mapping matrix initialized
nonzero_indices = nchoosek(1:K,N);
temp = eye(N,N);
flat = [];
for i=1:J
    V(nonzero_indices(i,:),:,i)= temp;
    %flat = [flat V(:,:,i)];
end
%%
ita = 0.005;
mod = 0.01;
%mod = 0;
result = [];
mutcoh = [];
num = 200; %number of simulation
count = 1;
%dicts = zeros(N,M*J,num); %cadidate dictionaries
for n=1:num
    for ii=1:length(ita)
        for jj=1:length(mod)
            %             init = randn(N,M*J)+1i*randn(N,M*J);
            init = randn(N,M)+1i*randn(N,M);
            
            B = sqrt(Es)*diag(1./vecnorm(init));
            init = init*B; %normalization
            H = Es*eye(M,M);
            iter1 = 10000;
            iter2 = 100;

            mu = sqrt((M-N)/(N*(M-1)));
            mu_mod = mu + mod(jj);
            
            for i1=1:iter1
                G = init'*init;
                for i=1:size(G,1)
                    for j=1:size(G,2)
                        
                        if i~=j
                            if abs(G(i,j)) > mu_mod*Es
                                H(i,j) = (G(i,j)/abs(G(i,j)))*mu_mod*Es;                                
                            else
                                H(i,j) = G(i,j);
                            end
                        end
                        
                    end
                    
                end
                abs(H);
                for i2=1:iter2
                    % column normalization
                    B = sqrt(Es)*diag(1./vecnorm(init));
                    init = init*B;
                    des = 4*init*init'*init-2*init*H-2*init*H';
               
                    init = init - ita(ii)*des;
                end
            end
            
            %vecnorm(init);
            
%             B = diag(1./vecnorm(init));
%             init = init*B;
            %result(ii,jj) = max(abs(GG(:)));
            GG = init'*init;
            GG = GG - diag(diag(GG));
            pdm = zeros(M,M);
            %PDM Calculation
            for icol=1:M
                for jcol=1:M
                    temp = prod(abs(init(:,icol)-init(:,jcol)));
                    if icol==jcol
                        pdm(icol,jcol) = inf;
                    else
                        pdm(icol,jcol) = temp;
                    end
                end
            end
            
            result = [result min(pdm(:))];
            mutcoh = [mutcoh max(abs(GG(:)))];
            dicts(:,:,n) = init;
            %             GG = GG - diag(diag(GG));
%             if min(abs(real(init(:))))>0.2 && min(abs(imag(init(:))))>0.2
%                 %result =[result max(abs(GG(:)))];
%                 result = [result min(pdm(:))];
%                 dicts(:,:,count) = init;
%                 count = count + 1;
%             end
        end
    end
end

%%
[val ind] = max(result)
[val1 ind1] = min(mutcoh)

non_sparse = dicts(:,:,ind1);
% B = diag(1./vecnorm(non_sparse));
% non_sparse = non_sparse*B; %normalization
final = [];
for i=1:J
    %     final = [final V(:,:,i)*non_sparse(:,(i-1)*M+1:i*M)];
    final = [final V(:,:,i)*non_sparse];
end

%save constellation_5SE.mat non_sparse
% save dict_mother.mat final