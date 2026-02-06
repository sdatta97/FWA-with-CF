function rate_dl = compute_link_rates_SISO(params,channel_dl)
M = params.numGNB;
K_FWA = params.numCPE;
K = params.numUE + params.numCPE;
BW = params.Band;
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
p_d = params.rho_tot; % 1*K;
D = params.D;

%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);
%Construc the above array and cells
for k = 1:K
    servingBSs = find(D(:,k)==1);
    NoservingBSs = find(D(:,k)==0);
    
    Serv{k} = servingBSs;
    NoServ{k} = NoservingBSs;
end

%% initialization of c
DS_Cell = zeros(K-K_FWA,1);
noise_Cell = abs(sqrt(0.5)*(randn(K-K_FWA,1) + 1j*randn(K-K_FWA,1))).^2;
rate_dl = zeros(K-K_FWA,1);
for k = 1:K-K_FWA
    DS_Cell(k) = p_d*(norm(channel_dl(Serv{k},k,:),'fro'))^2;
    rate_dl(k) = rate_dl(k) + (BW/(K-K_FWA))*TAU_FAC*log2(1+DS_Cell(k)/noise_Cell(k));
end
end