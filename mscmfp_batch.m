close all;
tic;
% main_modes;
% sources(1,:) = [1 1 1 3 3 3 3 5 5 5]/6;
% sources(2,:) = [5 9 13 3 7 11 15 5 9 13]/18;

%DEFAULTS:
clear; samplesize = 1e1 ; load states/state;

% pause(60*5);


% %%%%%%%%%%%%%%%%%% Attenuation Maps (1)
clear rhat;
% Mlist = 10;
nlevel = 0;

scenario = 'single';
sources = (2*eye(2)+1)/4;
al = [10; 1];
Mlist = 37;
clear G nr;
mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
state.map_single = reshape(normByCols(Up'*PGr)./normByCols(PGr),lzs,lrs);
state.rss = rs;
state.zs = zs;
% targets0(:,1) = targets(:,1);
state.targets_single = rhat;
erad0(:,1) = erad2;

clear G nr sources sources0 al;
scenario = 'coherent';
sources(1,:) = [-1 0 0 1]/3;
sources(2,:) = [0 -1 1 0]/3;
al = [2; 2; 2; 1];
sources = sources+1/2;
sources0 = sources;
nr = 6;
mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
state.map_coherent = reshape(normByCols(Up'*PGr)./normByCols(PGr),lzs,lrs);
state.rsc = rs;
state.targets_coherent = rhat;
erad0(:,2) = erad2;
save('data/attenuation','state','Mlist','erad0');
clear; samplesize = 1e3;  return;





% %%%%%%%%%%%%%%%%%% Covar Matrices (2 & 3)
Mlist = 37;
erqlist = (1:100)*.05;

clear G nr;
scenario = 'single';
mscmfp_simulate;

cS = zeros(M,100);
for k=1:length(erqlist)
    eradQk = erqlist(k);
    dx = sqrt(sum((DE*(grid_pts - rhat(:,k7)*ones(1,lzs*lrs))).^2));
    anx = find( dx < eradQk);
    Q = PG(:,anx)*PG(:,anx)';
    S = svd(Q);
    cS(:,k) = cumsum(S)/sum(S);
end;
cSe_single = cS;
Q = PG*PG';
S = svd(Q);
cSa_single = cumsum(S)/sum(S);

clear G nr ;
scenario = 'coherent';
mscmfp_simulate;
PG = permute(PG,[1 3 2]);
PG = reshape(PG,M*K,lrs*lzs);

cS = zeros(M*K,100);
for k=1:length(erqlist)
    eradQk = erqlist(k);
    dx = sqrt(sum((DE*(grid_pts - targets(:,1)*ones(1,lzs*lrs))).^2));
    anx = find( dx < eradQk);
    PGk = PG(:,anx);
    Q = PGk*PGk';
    S = svd(Q);
    cS(:,k) = cumsum(S)/sum(S);
    k,toc,pause(.01);
end;
cSe_coherent = cS;
Q = PG*PG';
S = svd(Q);
cSa_coherent = cumsum(S)/sum(S);

save('data/covar','cSe_single','cSa_single','cSe_coherent','cSa_coherent','erqlist');
clear; samplesize = 1e3;
return;






%%%%%%%%%%%%%%%%%%Region of Uncertainty (5)
sources = zeros(2);
sources(:,1) = [0.5; 10/180];
S0 = size(sources,2);
% num_iter = 10*S0; %number of OMP-style iterations
al = [10;1];
uk=[];

% scenariolist = {'single','coherent'};
% for k=1:length(scenariolist)
clear G nr;
clear schmall;
scenario = 'single';
nlevel=10^-1 * min(al)/norm(al);
erad2 = [36; 3];
Mlist = [37];
xx = 1:120;
yy = 1:90;
[XX,YY]=meshgrid(xx,yy);
xxyy = [XX(:)/120 YY(:)/180]'; %10800
nr = 3;
for k2=1:length(xxyy)
    sources(:,2) = xxyy(:,k2);
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    schmall(:,k2) = scm;
    title(num2str([k2/length(xxyy), toc/60]));
    pause(0.001);
end;
save(['data/RoU_' scenario],'rs','zs','Mlist','erad2','uk','schmall','xxyy','al','targets','xx','yy');
clear; samplesize = 1e3;

scenario = 'coherent';
al = [10; 1];
nlevel=10^(-0) / norm(al) * min(al);
erad2 = [12; 3];
% num_iter = 10;
% Mlist = [4 10];
Mlist = [37];
xx = 1:120;
yy = 1:60;
nr = 10;
% [XX,YY]=meshgrid(xx(20:5:end),yy(1:5:end));
[XX,YY]=meshgrid(xx,yy);
xxyy = [XX(:)/120 YY(:)/180]'; %10800
for k2=1:length(xxyy)
    sources(:,2) = xxyy(:,k2);
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    schmall(:,k2) = scm;
    rep = toc/60*[1, length(xxyy)/k2];
    title(num2str(rep));
    pause(0.001);
    k2
end;
save(['data/RoU_' scenario],'rs','zs','Mlist','erad2','schmall','xxyy','al','targets','xx','yy');
clear; samplesize = 1e3;

return;





% %%%%%%%%%%%%%%%%%% Single-Frequency(6)
clear G nr;
scenario = 'single';
Mlist = [2:37];
sources = (2*eye(2)-1)/3+1/2;
sources0 = sources;
scenario = 'single';
al = ones(2,1);
nlevel=10^-1 * min(al)/norm(al);
schmall = zeros(samplesize,length(Mlist));
reall = schmall;
k1 = 1;
for k2=1:samplesize
    sources = sources0 + randn(size(sources))/120;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    reall(k2,:,k1) = re; %residual energy
    schmall(k2,:,k1) = scm;
    title(num2str([k2, toc/60]));
    pause(0.001);
end;
squeeze(mean(schmall>1))

save(['data/multi_' scenario],'rs','schmall','zs','targets','Mlist','erad2','al','reall','nr');
clear; samplesize = 1e3; 






%%%%%%%%%%%%%%%%%% Coherent Final (7a)
scenario = 'coherent';
sources(1,:) = [-1 0 0 1]/3;
sources(2,:) = [0 -1 1 0]/3;
sources = sources+1/2;
sources = sources+randn(size(sources))/90;
Mlist = [4];
mscmfp_simulate;
save(['data/hex10_final_' scenario],'rs','zs','Mlist','scm','al','ht','targets','rhat');
return;

%%%%%%%%%%%%%%%%%% Coherent Tail (7b)
scenario = 'coherent';
sources(1,:) = [-1 0 0 1]/3;
sources(2,:) = [0 -1 1 0]/3;
sources = sources+1/2;
sources0 = sources;
Mlist = [4 37];
nr = 8;
schmall = zeros(samplesize,length(Mlist));
reall = schmall;
for k2=1:samplesize
    sources = sources0 + randn(size(sources))/120;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    title(num2str([k2, toc/60])); pause(.01);
    schmall(k2,:) = scm;
end;
% load(['data/hex10_' scenario]);
uk = sort(schmall,'descend');
R2.uk = uk;
save(['data/hex10_' scenario],'rs','zs','Mlist','erad2','uk','schmall');
clear; samplesize = 1e3;




%%%%%%%%%%%%%%%%%% ModErr (8)
Mlist = [4 37];
al = [1; 1]; %track each one individually
scenario = 'coherent';
sources0 = (2*eye(2)+1)/4;
schmall = zeros(samplesize,length(Mlist),6);
erad2 = [1;1];
for k3=1:6
    eval(['load states/state_141_160_' num2str(1518+2*k3) '.mat;']);
    for k2=1:samplesize
        sources = sources0+randn(2)/120;
        mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
        schmall(k2,:,k3) = scm;
    end;
    clear G;
    k3,toc,pause(1);
end;
uk = sort(schmall,'descend');
R2.uk = uk;
save(['data/moderr_' scenario],'rs','zs','Mlist','erad2','uk','schmall');
clear; samplesize = 1e3; 




%%%%%%%%%%%%%%%%%% Noise tests (9)
scenario = 'single';
sources0 = eye(2)/2+1/4;
Mlist = [37];
al = ones(2,1);
schmall = zeros(samplesize,5);
for k1=1:5
dBlevel(k1) = 0.1*10^(-(k1-3)/4);
nlevel = dBlevel(k1) * min(al)/norm(al);
for k2=1:samplesize               
    sources = sources0 + randn(size(sources0))/120;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    title(num2str([k2, toc/60])); pause(.01);
    schmall(k2,k1) = scm;
end;
end;
save(['data/noise_' scenario],'Mlist','schmall','dBlevel');
clear G nr;
clear; samplesize = 1e3;

scenario = 'coherent';
sources(1,:) = [-1 0 0 1]/3;
sources(2,:) = [0 -1 1 0]/3;
sources = sources+1/2;
sources0 = sources;
Mlist = [37];
al = ones(4,1);
schmall = zeros(samplesize,5);
for k1=1:5
% dBlevel(k1) = 10^(-k1/5);
dBlevel(k1) = 10^(-(k1-3)/4);
nlevel = dBlevel(k1) * min(al)/norm(al);
for k2=1:samplesize               
    sources = sources0 + randn(size(sources))/120;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    title(num2str([k2, toc/60])); pause(.01);
    schmall(k2,k1) = scm;
end;
end;
save(['data/noise_' scenario],'Mlist','schmall','dBlevel');
clear G nr;
clear; samplesize = 1e3;



% %%%%%%%%%%%%%%%%%% Point Nulling (10)
scenario = 'coherent';
Mlist = [4 37];
sources(1,:) = [-1 0 0 1]/12;
sources(2,:) = [0 -1 1 0]/12;
sources = sources+1/2;
% sources = rand(2,4);
sources0 = sources;
schmall = zeros(samplesize,length(Mlist));
reall = schmall;
al = ones(10,1);
nlevel=10^-0 * min(al)/norm(al);
for k2=1:samplesize
    sources = sources0+rand(size(sources0))/120;
    
    nr = 3;
    eradQ=  .01;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    schmall(k2,1:2) = scm;
    reall(k2,1:2) = re;
    
    nr = 6;
    eradQ = 0.4;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    schmall(k2,3:4) = scm;
    reall(k2,3:4) = re;
    disp(num2str(median(schmall(1:k2,:))));
    
    disp(' ');
    title(num2str([k2, toc/60])); pause(.01);
end;
save(['data/point_nulling_' scenario],'rs','zs','Mlist','schmall','reall');
sort(schmall)
return;
clear G nr;
clear; samplesize = 1e3; load states/state;



%%%%%%%%%%%%%%%%%% Clean (11)
eradQ=  0.01;
flags = '';
scenario = 'coherent';
Mlist = [4 37]; 
sources(1,:) = [-1 0 0 1]/3;
sources(2,:) = [0 -1 1 0]/3;
sources = sources+1/2;
sources0 = sources;
schmall = zeros(samplesize,length(Mlist)*2);
reall = schmall;
for k2=1:samplesize
    sources = sources0+rand(size(sources0))/120;
    
    nr = 8;
    eradQ = 1/1;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    schmall(k2,3:4) = scm;
    reall(k2,3:4) = re;

    lgain = 0.9;
    mscmfp_clean;  %%%%%% ACTION HERE %%%%%
    nx=1:2;
    schmall(k2,nx) = scm;

    disp(num2str(median(schmall(1:k2,:))));
    title(num2str([k2, toc/60])); pause(.01);
end;
save(['data/clean_' scenario],'rs','zs','Mlist','schmall','reall');
clear G nr;
clear; samplesize = 1e3;


% %%%%%%%%%%%%%%%%%% Global Search (12)
flags = 'global';
scenario = 'coherent';
Mlist = [4];
sources0 = eye(2)/2+1/4;
schmall = zeros(samplesize,1);
schmallg = schmall; reall = schmall; regall = schmall;
nr = 8;
for k2=1:samplesize
    sources = sources0+rand(2)/120;
    mscmfp_simulate;  %%%%%% ACTION HERE %%%%%
    schmall(k2,1)  = scm;
    schmall(k2,2) = scmg;
    reall(k2,1) = re;
    reall(k2,2) = reg;
    regall(k2,:) = reg;
    title(num2str([k2, toc/60])); pause(.01);
    disp(num2str(median(schmall(1:k2,:))));
end;
save(['data/global_' scenario],'rs','zs','Mlist','schmall','reall');
%save
clear flags G;
clear; samplesize = 1e3;



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%NETHER REGION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% %%%%%%%%%%%%%%%%%% Lipchitz Profile (12)
scenario = 'coherent';
mscmfp_simulate;
G = reshape(G,N*K,lrs*lzs);
G = G./(ones(N*K,1)*normByCols(G));
mscmfp_lipchitz;  %%%%%% ACTION HERE %%%%%
save(['data/lipchitz_' scenario],'minm','maxm','maxt','maxt2');
clear G nr;
pause;

scenario = 'single';
mscmfp_simulate;
G = squeeze(G(:,10,:));
G = G./(ones(N,1)*normByCols(G));
mscmfp_lipchitz;  %%%%%% ACTION HERE %%%%%
save(['data/lipchitz_' scenario],'minm','maxm','maxt','maxt2');
clear G nr;
clear; samplesize = 1e1;  return;



