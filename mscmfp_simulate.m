%all (x,y) coordinates are in the format: range/depth

% Defaults
if ~exist('scenario');    scenario='coherent';                 end;
if ~exist('phi_type');    phi_type='QR-Gaussian';              end;
if ~exist('Mlist');       Mlist = 37;                          end;
if ~exist('sources');     sources = (2*eye(2)+1)/4;            end;
if ~exist('numSrcs');     numSrcs = size(sources,2);           end;
if ~exist('num_iter');    num_iter = numSrcs*10;               end;
if ~exist('al');          al = ones(numSrcs,1);                end;
if ~exist('nlevel');      nlevel=10^-0 * min(al)/norm(al);     end;
if ~exist('eradQ');       eradQ = 1;                           end;
if ~exist('sing_factor'); sing_factor = 0;                     end;
if ~exist('flags');       flags = '';                          end;
if ~exist('erad2');       erad2 = [72; 6];                     end;

 %[36;3];

% if ~exist('psi)');        load states/state_141_160_1520.mat;  end;

if strcmp(scenario,'coherent'),
if ~exist('nr');          nr = numSrcs*2;                      end;
end;

numSrcs = size(sources,2);
numSrcs_os = numSrcs+0;


rs = 5000+(1:(240))*12;  %5000+(1:(120))*6;
zs = 10+(1:(90))*2;             %Source Depth (m)  (10+(1:(180))*1)
zr = 10:5:190;                  %Receiver Depth (m)
dls = [mean(diff(rs)); mean(diff(zs))];
DE = diag(erad2.^(-1));

lrs = length(rs);
lzs = length(zs);
lzr = length(zr);
lf = length(freq);
if strcmp(scenario,'single'),      flist = lf/2;               end;
if strcmp(scenario,'incoherent'),  flist = 1:lf;               end;
if strcmp(scenario,'coherent'),    flist = 1:lf;               end;

cmap = colormap('gray');
colormap(flipud(cmap));

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
    if ~strcmp(scenario,'coherent') %set explicitly for coherent
        nr = min(6,round(M/2));
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
    
    if findstr(flags,'global')
        mscmfp_global;
        [yv,dx] = minimax_assignment(DE*targets,DE*rhat);
        scmg(ii) = yv;
        regg(ii) = reg;
        
        %Status update / display
        disp(['With M=' num2str(M) ' tests, scm(M) = ' num2str(scmg(ii)), ', ' num2str(toc/60) ' minutes have elapsed. (global)']);
        pause(.01);
    end;

    %Round-robin greedy localization
    beta = zeros(numSrcs_os,KF);
    for k5=1:num_iter
        eradQk = eradQ;% * .1^(k5/num_iter);
        Up = [];
        ks = mod(k5-1,numSrcs_os)+1;
        rhat(:,ks) = 0;
        PGhat(:,:,ks) = 0;
        
        %Generate the covariance matrix of compressed Green's functions
        R  = zeros(M,M,KF);
        Rc = zeros(M*KF,M*KF);
        for k7 = 1:numSrcs_os
            if (norm(rhat(:,k7))~=0)
            dx = sqrt(sum((DE*(grid_pts - rhat(:,k7)*ones(1,lzs*lrs))).^2));
            anx = find( dx < eradQk);
            Gl = zeros(M,length(anx),KF);
            for kf=1:KF
                k = flist(kf);
                temp = PG(:,anx,kf);
                Gl(:,:,kf) = temp;
                R(:,:,kf) = R(:,:,kf) + temp*temp';
            end;
            Gl = reshape(permute(Gl,[1 3 2]),M*KF,length(anx));
            Rc = Rc + Gl*Gl';
            end; %if
        end; %k7
        
        %Construct the compressed ambiguity function
        if strcmp(scenario,'coherent')
            den2 = sum(nPG2);
            UUPY = PY;
            if norm(rhat,'fro')~=0
                if M*KF<200
                    [U,S,V] = svd(Rc);
                else %fast approximation to U
                    [U,S,V] = svd(Rc*randn(M*KF,max(100,nr*2)));
                end;
                Up = U(:,(nr+1):end);
                UUPY = Up*Up'*PY(:); %creative destruction
                den2 = den2 - sing_factor*normByCols(U(:,1:nr)'*PGr).^2; %*C3*
                %re-introduce this orthogonalizing line
            end;
            ht = abs(UUPY(:)'*PGr).^2 ./ den2 / norm(PY(:))^2;
        else
            ht = 0;
            for kf=1:KF
                if norm(rhat,'fro')==0
                    Up = eye(M);
                else
                    [U,S,V] = svd(R(:,:,kf));
                    Up = U(:,(nr+1):end);
                end;
                PGk = Up'*PG(:,:,kf);
                PYk = Up'*PY(:,kf);
                %                 ht = ht + abs(PYk'*PGk).^2./nPG2(kf,:);
                ht = ht + abs(PYk'*PGk).^2./normByCols(PGk).^2;
%                 ht = ht + abs(PYk'*PGk).^2./normByCols(PG(:,:,kf)).^2;
            end;
            ht = ht / norm(PY(:))^2;
        end;
        ht = abs(ht);
        imagesc(rs,zs,reshape(ht,lzs,lrs));
        hold on;
        plot(targets(1,:),targets(2,:),'ro','MarkerSize',12);
        plot(rhat(1,:),rhat(2,:),'kx','MarkerSize',12);
%         plot(rhat_prev(1,:),rhat_prev(2,:),'c+','MarkerSize',12);
        hold off; pause(.01);
%         cmap = colormap(gray);
%         colormap(flipud(cmap));
%         xlabel('Range (m)');
%         ylabel('Depth (m)');
%         colorbar;
        
        %drumroll....
        [rek1,cmnx] = max(ht);
        rhat(:,ks) = grid_pts(:,cmnx);
        PGhat(:,:,ks) = squeeze(PG(:,cmnx,:));
        if findstr(flags,'movie')
            frame= getframe;
            mov = addframe(mov,frame.cdata(1:336,1:432,:));
        end;

        
        %Current beta estimate
        %the only take-away here is the residual energy... is this necessary?
        if strcmp(scenario,'coherent')
            A = reshape(PGhat,M*KF,numSrcs_os);
            b = reshape(PY,M*KF,1);
            beta = pinv(A)*b;
            re(ii) = norm(A*beta-b)^2 / norm(b)^2;
            beta = beta*ones(1,KF);
        else
            re(ii) = 0;
            for kf=1:KF
                A = squeeze(PGhat(:,kf,:));
                b = PY(:,kf);
                beta(:,kf) = pinv(A)*b;
                re(ii) = re(ii) + norm(A*beta-b)^2 / norm(b)^2 / KF;
            end;
        end;
        if (ks==numSrcs_os), rhat_prev = rhat; end;
    end; %k5
    [y,nx] = sort(normByCols(beta'),'descend');
    rhat = rhat(:,nx(1:numSrcs));
    [yv,dx] = minimax_assignment(DE*targets,DE*rhat);
    scm(ii) = yv;
    
    %Status update / display
    disp(['With M=' num2str(M) ' tests, scm(M) = ' num2str(scm(ii)), ', ' num2str(toc/60) ' minutes have elapsed. Residual energy: ' num2str(re(end))]);
    pause(.01);
end; %ii
%%%%%%%%%%%%% -------- SIMS FINISHED ---------- %%%%%%%%%%

% %SPARE HOURS
% hh = mod(now*24,24);
% if ( ((hh>9)&(hh<11)) | (hh>21) )
%     pause(300); %run this on <10% duty cycle during peak usage
% end;

return;
