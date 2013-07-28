%all (x,y) coordinates are in the format: range/depth

% Defaults
if ~exist('scenario');    scenario='single';                   end;
if ~exist('phi_type');    phi_type='QR-Gaussian';              end;
if ~exist('Mlist');       Mlist = [6];                         end;
if ~exist('numSrcs');     numSrcs = size(sources,2);           end;
if ~exist('al');          al = ones(numSrcs,1);                end;
if ~exist('nlevel');      nlevel=10^-1 * min(al)/norm(al);     end;
if ~exist('eradQ');       eradQ = 1/1;                         end;
if ~exist('flags');       flags = '';                          end;
if ~exist('lgain');       lgain = 0.9;                         end;


numSrcs = size(sources,2);
numSrcs_os = numSrcs+0;

erad2 = [72; 6]; %[36;3];
rs = 5000+(1:(240))*12;  %5000+(1:(120))*6;
if strcmp(scenario,'coherent'),    erad2 = [24;6];             end; %[12;3];
if strcmp(scenario,'coherent'),    rs = 5000+(1:(240))*12;     end; %5000+(1:(120))*2;

zs = 10+(1:(90))*2;             %Source Depth (m)  (10+(1:(180))*1)
zr=10:5:190;                    %Receiver Depth (m)
dls = [mean(diff(rs)); mean(diff(zs))];
DE = diag(erad2.^(-1));

lrs = length(rs);
lzs = length(zs);
lzr = length(zr);
lf = length(freq);
if strcmp(scenario,'single'),      flist = lf/2;               end;
if strcmp(scenario,'incoherent'),  flist = 1:lf;               end;
if strcmp(scenario,'coherent'),    flist = 1:lf;               end;

%aliases
N = lzr; 
K = length(freq);
KF = length(flist);

% Indexed by kf: PY, PG, PGhat, rhat
% Indexed by k:  Yom, G, rGg, Ggt, 

trs = (0*zs'+1)*rs;
tzs = zs'*(0*rs+1);
grid_pts = [trs(:) tzs(:)]';

%constrain to the unit square
sources(sources< 0) =  0;
sources(sources> 1) =  1;
SS = [mean(diff(rs))*length(rs); mean(diff(zs))*length(zs)];
targets = diag(SS)*(sources) + [rs(1)-dls(1); zs(1)-dls(2)]*ones(1,numSrcs);

if ~exist('G') %this will get me in trouble later, but should be fine for now.
    G = greensG_mode(psi,z,N_modes,modes,rho_w,zr,grid_pts);
    G = permute(G,[1 3 2]) / max(abs(G(:)));
end;

%G is only computed once, there is opportunity to reload the state
%information at every stage here. It is possible to then run the modeling
%error simulation here.

%reload for the modeling error sims, unnecessary in general
load states/state_141_160_1520.mat; 

Ggt = zeros(N,K,numSrcs);
for l=1:numSrcs
    temp = greens_mode(psi,z,N_modes,modes,rho_w,zr,targets(:,l));
    temp2 = reshape(temp,K,N).';
    Ggt(:,:,l) = temp2/norm(temp2(:,round(lf/2))); %normalization placeholder
end;
Ggt = permute(Ggt,[1 3 2]);
rGg = zeros(N,K);
for k=1:K
%     rGg(:,k) = Ggt(:,:,k)*(al.*exp(j*2*pi*rand(size(al))));
    rGg(:,k) = Ggt(:,:,k)*(al);
end;
if strcmp(scenario,'coherent')
    temp = reshape(permute(Ggt,[1 3 2]),K*N,numSrcs);
    ipp = abs(temp(:,1)'*temp(:,2))/norm(temp(:,1))/norm(temp(:,2));
else;
    ipp = abs(Ggt(:,1,10)'*Ggt(:,2,10))/norm(Ggt(:,1,10))/norm(Ggt(:,2,10)); %dirty
end;

num  = 0;
num2 = 0;
den2 = 0;

nn = (randn(N,K)+j*randn(N,K))/sqrt(2*N*K);
Yom = rGg + nn * nlevel*norm(rGg(:));

%%%%%%%%%%%%% -------- COMPRESSED MFP ---------- %%%%%%%%%%
scm = zeros(1,length(Mlist)); %initialize output
for ii=1:length(Mlist)
    M = Mlist(ii);  %number of test-measurements
    if strcmp(scenario,'coherent')
        nr = min(round(1.5*numSrcs_os),20);
    else
        nr = min(floor(M/2),5);
    end;
    loclist = zeros(1,1);
    rhat = zeros(2,numSrcs_os);
    rhat_prev = rhat;
    
    %Generate Phi
    Phi = zeros(M,N,KF);
    if findstr(phi_type,'Bernoulli'),     Phif = sign(randn(lzr,lzr,KF));
    elseif findstr(phi_type,'Gaussian'),  Phif =      randn(lzr,lzr,KF);
    end;
    PG = zeros(M,lzs*lrs,KF);
    PY = zeros(M,KF);
    nPG2 = zeros(KF,lrs*lzs);
    for kf=1:KF
        k = flist(kf);
        Phik = Phif(:,:,kf);
        if findstr(phi_type,'QR'),
            [Phik,R]=qr(Phik);
        end;
        Phi(:,:,kf)=Phik(:,1:M)';
        PY(:,kf)   = Phi(:,:,kf)*Yom(:,k);
        PG(:,:,kf) = Phi(:,:,kf)*squeeze(G(:,k,:));
        nPG2(kf,:) = sum(abs(PG(:,:,kf)).^2);
    end;
    PGr = reshape(permute(PG,[1 3 2]),M*KF,lrs*lzs);
    PGhat = zeros(M,KF,numSrcs_os);
    
    %Round-robin greedy localization
    beta = zeros(numSrcs_os,KF);
%     for k5=1:num_iter
    for k5=1:numSrcs
        eradQk = eradQ;% * .1^(k5/num_iter);
        Up = [];
        ks = mod(k5-1,numSrcs_os)+1;
        rhat(:,ks) = 0;
        PGhat(:,:,ks) = 0;
        
        %Construct the compressed ambiguity function
        if strcmp(scenario,'coherent')
            den2 = sum(nPG2);
            ht = abs(PY(:)'*PGr).^2 ./ den2 / norm(PY(:))^2;
        else
            ht = 0;
            for kf=1:KF
                PGk = PG(:,:,kf);
                PYk = PY(:,kf);
                ht = ht + abs(PYk'*PGk).^2./normByCols(PGk).^2;
            end;
            ht = ht / norm(PY(:))^2;
        end;
        ht = abs(ht);
        imagesc(rs,zs,reshape(ht,180,120));
        hold on;
        plot(targets(1,:),targets(2,:),'ro','MarkerSize',12);
        plot(rhat(1,:),rhat(2,:),'wx','MarkerSize',12);
        hold off; pause(.01);
        pause(.01);
        
        %drumroll....
        [rek1,cmnx] = max(ht);
        rhat(:,ks) = grid_pts(:,cmnx);
        PGhat(:,:,ks) = squeeze(PG(:,cmnx,:));
        
        dx = sqrt(sum((DE*(grid_pts - rhat(:,ks)*ones(1,lzs*lrs))).^2));
        anx = find( dx < eradQk);
        Gl = zeros(M,length(anx),KF);
        
        if strcmp(scenario,'coherent')
            for kf=1:KF
                k = flist(kf);
                temp = PG(:,anx,kf);
                Gl(:,:,kf) = temp;
            end;
            Gl = reshape(permute(Gl,[1 3 2]),M*KF,length(anx));
            Rc = Gl*Gl';
            Rc = Rc / trace(Rc);
            
            PGk = PGr(:,cmnx);
            PGk = PGk/norm(PGk);
            PY = PY(:) - Rc*PY(:);
%             PY = PY(:) - PGk*PGk'*PY(:);
        else
            PGk = PGhat(:,:,ks);
            for kf=1:KF
                k = flist(kf);
                temp = PG(:,anx,kf);
                R = temp*temp';
                R = R/trace(R);
                
                PGk = PG(:,cmnx,kf);
                PGk = PGk/norm(PGk);
                PY(:,kf) = PY(:,kf) - R*PY(:,kf);
            end;
        end;
    end; %k5
    [yv,dx] = minimax_assignment(DE*targets,DE*rhat);
    scm(ii) = yv;
    
    %Status update / display
    disp(['With M=' num2str(M) ' tests, scm(M) = ' num2str(scm(ii)), ', ' num2str(toc/60) ' minutes have elapsed.']);
    pause(.01);
end; %ii
%%%%%%%%%%%%% -------- SIMS FINISHED ---------- %%%%%%%%%%


%SPARE HOURS
hh = mod(now*24,24);
if ( ((hh>9)&(hh<11)) | (hh>21) )
%     pause(60); %run this on 10% duty cycle during peak usage
end;

return;

% if ((mod(now,1)*24>10) & (mod(now,1)*24<21) ), pause(3600), end;