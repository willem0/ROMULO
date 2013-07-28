% cmfp_figGen.m
clear;
tic;
close all;

set(0,'DefaultAxesFontsize',18);
cmap = colormap('hot'); cmap = flipud(cmap);

load states/state_141_160_1520.mat
phi_type = 'QR-Gaussian';

shg;


% nlevel= 0.1 * min(al)/norm(al);






%%%%%%%%%%%%%%%%%% Point-nulling (9)
scenario = 'coherent';
load(['data/clean_' scenario]);
samplesize = length(schmall);
pp = linspace(1,1e-3,samplesize);
u = schmall(:,2);
scenario = 'coherent';
load(['data/point_nulling_' scenario]);
uk = sort([schmall u]);
semilogy(uk(:,1),pp,'--',uk(:,3),pp,'-',uk(:,2),pp,'--',uk(:,4),pp,'-');
xlabel('distance d');
ylabel('P_M\{error > d\}');
axis([0 1 1e-2 1]);
legend({'M=4 (point)','M=4 (robust)','M=37 (point)','M=37 (robust)'});
clipboard('copy',['point_nulling_' scenario]);
pause;;;



%%%%%%%%%%%%%%%%%% Region of Uncertainty (5)

scenario = 'single';
load(['data/RoU_' scenario]);
% xx = unique(xxyy(1,:));
% yy = unique(xxyy(2,:));
UR = reshape(schmall(1,:),length(yy),length(xx));
UR2 = mscmfp_sigmoid((UR-1)*5);
UR3 = [mscmfp_gsmooth(mscmfp_gsmooth(UR2',1)',1); zeros(90,120)];
theta = linspace(0,2*pi);
imagesc(5000+24*(1:120),10+(1:90),1-UR3);
hold on;
plot(targets(1,1)+cos(theta)*erad2(1),targets(2,1)+sin(theta)*erad2(2),'g');
% plot(5360+cos(theta)*erad2(1),11+sin(theta)*erad2(2),'w');
hold off;
colormap(flipud(cmap));
xlabel('range(m)');
ylabel('depth (m)');
clipboard('copy',['RoU_' scenario]);
colorbar;
pause;

scenario = 'coherent';
load(['data/RoU_' scenario]);
% xx = unique(xxyy(1,:));
% yy = unique(xxyy(2,:));
UR = reshape(schmall(1,:),length(yy),length(xx));
UR2 = mscmfp_sigmoid((UR-1)*5);
UR3 = [mscmfp_gsmooth(mscmfp_gsmooth(UR2',1)',1); zeros(120)];
imagesc(5000+24*(1:120),10+(1:180),1-UR3);
hold on;
% erad2 = erad2*10;
plot(targets(1,1)+cos(theta)*erad2(1),targets(2,1)+sin(theta)*erad2(2),'g');
hold off;
colorbar;
xlabel('range(m)');
ylabel('depth (m)');
colormap(flipud(cmap));
clipboard('copy',['RoU_' scenario]);

pause;;;





% %%%%%%%%%%%%%%%%%% Modeling Error (13)
scenario = 'coherent';
load(['data/moderr_' scenario]);
% N = length(scml);
% for k=1:N
%     mfp(k)  = (mean(scml(k).snm));
%     cmfp(k) = (mean(scml(k).scm));
% end;
mfp_cmfp = squeeze(median(schmall(1:1000,:,:)));
plot(0:2:10,mfp_cmfp(2,:),'-+',0:2:10,mfp_cmfp(1,:),'-o');
xlabel('modeling error (m/s)');
ylabel('localization error (m)');
legend({'MFP','cMFP'});
clipboard('copy',['moderr_' scenario]);
pause;
clear mfp cmfp;

scenario = 'single';
return;




%%%%%%%%%%%%%%%%%% Single Pair (6)
scenario = 'single';
% load(['data/multi_' scenario '_allM']);
load(['data/multi_' scenario]);
samplesize = length(schmall);
pp = linspace(1,1e-3,samplesize);
uk = sort(schmall);
ax(1) = find(Mlist==8);
ax(2) = find(Mlist==16);
ax(3) = find(Mlist==37);
% semilogy(uk(:,ax,2),pp);
semilogy(uk(:,ax(1),1),pp,'-',...
         uk(:,ax(2),1),pp,'-.',...
         uk(:,ax(3),1),pp,'--');
axis([0 1 1e-2 1]);
xlabel('distance d');
ylabel('P_M\{error > d\}');
for k=1:length(ax)
    legstr{k} = ['M=' num2str(Mlist(ax(k)))];
end;
legend(legstr);
clipboard('copy',['pair_' scenario]);
pause;

uk = mean(schmall(:,:,1)>1)';
% semilogy(Mlist,uk(:,1),'--s',Mlist,uk(:,2),'--^','MarkerSize',8);
semilogy(Mlist,uk);
axis([2 37 1e-2 1]);
xlabel('M');
ylabel('P_M\{error > d=1\}');
clipboard('copy',['pair_allM_' scenario]);
pause;;;







%%%%%%%%%%%%%%%%%% Global / Residual (11)
scenario = 'coherent';
load(['data/global_' scenario]);
samplesize = length(schmall);
pp = linspace(1,1e-3,samplesize);
uk = sort(schmall);
semilogy(uk(:,1),pp,uk(:,2),pp,'--');
xlabel('distance d');
ylabel('P_2\{error > d\}');
axis([0 1 1e-2 1]);
legend({'greedy','global'});
clipboard('copy',['global_' scenario]);
pause;

uk = sort(reall);
semilogy(uk(:,1),pp,uk(:,2),pp,'--');
xlabel('error e');
ylabel('P\{residual error (squared) > e\}');
axis([0 1 1e-2 1]);
legend({'greedy','global'});
clipboard('copy',['residual_' scenario]);
pause;;;





%%%%%%%%%%%%%%%%%% Clean (10)
scenario = 'coherent';
load(['data/clean_' scenario]);
samplesize = length(schmall);
pp = linspace(1,1e-3,samplesize);
uk = sort(schmall);
semilogy(uk(:,1),pp,'--',uk(:,3),pp,'-',uk(:,2),pp,'--',uk(:,4),pp,'-');
xlabel('distance d');
ylabel('P_M\{error > d\}');
axis([0 1 1e-2 1]);
legend({'M=4 (clean)','M=4','M=37 (clean)','M=37'});
clipboard('copy',['clean_' scenario]);
pause;;;




%%%%%%%%%%%%%%%%%%Intermediate Frames (4)
% scenario = 'single';
% sources(:,1) = [1;1]/4;
% sources(:,2) = [3;3]/4;
% erad2 = [36; 3];
% Mlist = [10];

scenario = 'coherent';
Mlist = [37];
al = ones(4,1);
nlevel= 1.0 * min(al)/norm(al);
sources(1,:) = [-1 0 0 1]/3;
sources(2,:) = [0 -1 1 0]/3;
sources = sources+1/2 + randn(size(sources))/90;
mscmfp_simulate;
xlabel('range(m)');
ylabel('depth (m)');
clipboard('copy',['multi_' scenario '_stage']);
pause;;;



%%%%%%%%%%%%%%%%%% Attenuation Map (1)
load('data/attenuation');
tt = linspace(0,2*pi);
erad0=[36 12; 3 3];


imagesc(state.rss,state.zs,(state.map_single));
xlabel('range(m)');
ylabel('depth (m)');
clipboard('copy','attenuation_single');
colormap(flipud(cmap));
colorbar;
hold on;
plot(cos(tt)*erad0(1,1)+state.targets_single(1,1),sin(tt)*erad0(2,1)+state.targets_single(2,1),'w');
hold off;
pause;

imagesc(state.rsc,state.zs,(state.map_coherent));
xlabel('range(m)');
ylabel('depth (m)');
clipboard('copy','attenuation_coherent');
colormap(flipud(cmap));
colorbar;
hold on;
for k=1:(length(state.targets_coherent))
plot(cos(tt)*erad0(1,2)+state.targets_coherent(1,k),sin(tt)*erad0(2,2)+state.targets_coherent(2,k),'w');
pause(1);
end;
hold off;
pause;;;



% %%%%%%%%%%%%%%%%%% Covar / Sigma Capture (2/3)

load('data/covar');

scenario = 'single';
semilogy(erqlist,1-cSe_single(1,:)','-',...
         erqlist,1-cSe_single(3,:)','-.',...
         erqlist,1-cSe_single(5,:)','--');
% semilogy(erqlist,1-cSe_single([1 3 5],:)');
axis([0 5 1e-4 1]);
xlabel('scale factor  \gamma');
ylabel('attenuation factor');
legend({'n_r=1','n_r=3','n_r=5'});
clipboard('copy',['attenuation_gamma_' scenario]);
pause;

% sigma_capture_single
xx = 0:length(cSa_single);
cSa_single = [0; cSa_single];
plot(xx,1-cSa_single,xx,1-xx*cSa_single(2),'--');
axis([0 15 0 1])
xlabel('n_r');
ylabel('attenuation factor');
legend({'attenuation','approximation'});
clipboard('copy',['sigma_capture_' scenario]);
pause;

scenario = 'coherent';
% semilogy(erqlist,1-cSe_coherent([3 6 9 12],:)');
semilogy(erqlist,1-cSe_coherent(3,:)','-',...
         erqlist,1-cSe_coherent(6,:)','-.',...
         erqlist,1-cSe_coherent(9,:)','--',...
         erqlist,1-cSe_coherent(12,:)','.-');
xlabel('scale factor  \gamma');
axis([0 5 1e-4 1]);
legend({'n_r=3','n_r=6','n_r=9','n_r=12'});
clipboard('copy',['attenuation_gamma_' scenario]);
pause;

% sigma_capture_coherent
xx = 0:length(cSa_coherent);
cSa_coherent = [0; cSa_coherent];
plot(xx,1-cSa_coherent,xx,1-xx*cSa_coherent(2),'--');
axis([0 2/cSa_coherent(2) 0 1])
xlabel('n_r');
legend({'attenuation','approximation'});
clipboard('copy',['sigma_capture_' scenario]);
pause;




%%%%%%%%%%%%%%%%%% Pattern (7a)
scenario = 'coherent';
load(['data/hex10_final_' scenario]);
imagesc(rs,zs,reshape(ht,180,120));
hold on;
plot(targets(1,:),targets(2,:),'ro','MarkerSize',12);
plot(rhat(1,:),rhat(2,:),'mx','MarkerSize',12);
hold off;
cmap = colormap(gray);
colormap(flipud(cmap));
xlabel('range (m)');
ylabel('depth (m)');
colorbar;
clipboard('copy',['hex10_final_' scenario]);
pause;






%%%%%%%%%%%%%%%%%% Noise (8)
scenario = 'single';
load(['data/noise_' scenario]);
samplesize = length(schmall);
pp = linspace(1,1e-3,samplesize);
uk = sort(schmall);
% semilogy(uk(:,1:end-1),pp);
semilogy(uk(:,1),pp,'-',...
         uk(:,2),pp,'-.',...
         uk(:,3),pp,'--',...
         uk(:,4),pp,'.-');
axis([0 1 1e-2 1]);
xlabel('distance d');
ystr = ['P_{' num2str(Mlist) '}\{error > d\}'];
ylabel(ystr);
for k=1:length(dBlevel)-1
    dBk = -db(dBlevel(k));
    legstr{k} = [num2str(dBk) ' dB'];
end;
legend(legstr);
clipboard('copy',['noise_' scenario]);
pause;

scenario = 'coherent';
load(['data/noise_' scenario]);
uk = sort(schmall);
semilogy(uk(:,1),pp,'-',...
         uk(:,2),pp,'-.',...
         uk(:,3),pp,'--',...
         uk(:,4),pp,'.-');
axis([0 1 1e-2 1]);
xlabel('distance d');
ystr = ['P_{' num2str(Mlist) '}\{error > d\}'];
ylabel(ystr);
for k=1:length(dBlevel)-1
    dBk = -db(dBlevel(k));
    legstr{k} = [num2str(dBk) ' dB'];
end;
legend(legstr);
clipboard('copy',['noise_' scenario]);
pause;;;













%%%%%%%%%%%%%%%%%% Coherent (7b)
scenario = 'coherent';
load(['data/hex10_' scenario]);
samplesize = length(schmall);
pp = linspace(1,1e-3,samplesize);

% uk = R2.uk;
uk = sort(schmall);
semilogy(uk(:,1),pp,'-',...
         uk(:,2),pp,'--');
axis([0 1 1e-2 1]);
xlabel('distance d');
ylabel('P_M\{error > d\}');

% for k=1:length(Mlist)
for k=1:length(Mlist)
    legstr{k} = ['cMFP: M=' num2str(Mlist(k))];
end;
legend(legstr);
clipboard('copy',['hex10_' scenario]);
pause;;;








%%%%%%%%%%%%%%%%%% Lipchitz (12)
scenario = 'coherent';
load(['data/lipchitz_' scenario]);
plot(maxt2,minm,'--',maxt2,maxm,'-.');
axis([0 30 1e-2 1]);
xlabel('physical distance (squared)');
ylabel('Green''s function distance (squared)');
clipboard('copy',['lipchitz_' scenario]);
legend({'lower bound','upper bound'});

scenario = 'single';
load(['data/lipchitz_' scenario]);
plot(maxt2,minm,'--',maxt2,maxm,'-.');
axis([0 30 1e-2 1]);
xlabel('physical distance (squared)');
ylabel('Green''s function distance (squared)');
clipboard('copy',['lipchitz_' scenario]);
legend({'upper bound','lower bound'});
pause;;;





%%%%%%%%%%%%%%%%%% CoSaMP Comparison
load data/CoSaMP_coherent.mat
N = 1000;
semilogy([1 1],[eps 1],'--', (R2.uk_2s),(1:N)/N, (R2.uk_1s),(1:N)/N, R2.uk_reprojection,(1:N)/N, R2.uk_null,(1:N)/N);
% semilogy([1 1],[eps 1],'--', (R2.uk_2s),(1:N)/N, (R2.uk_1s),(1:N)/N, R2.uk_reprojection,(1:N)/N, R2.uk_null,(1:N)/N, uk,(1:100)/100);
xlabel('distance d');
ylabel('P_4\{error > d\}');
legend({'target ellipse','2s','1s','reprojection','null','ms-cMFP'});
axis([0 4 .01 1])
return;



%%%%%%%%%%%%%%%%%%%%%%%% Sigma Capture (Obsolete)
scenario = 'single';
load(['data/sigcapture_' scenario],'S','Sc','Mlist');
xx = 1:length(S);
xc = 1:length(Sc);
% plot(xx,1-cumsum(S),xc,1-cumsum(Sc));
plot(xx,1-cumsum(S),xx,1-xx*S(1),'--');
xlabel('n_r');
ylabel('mean-square attenuation');
legend({'attenuation','approximation'});
axis([1 15 0 1])
clipboard('copy',['sigma_capture_' scenario]);
pause;

scenario = 'coherent';
load(['data/sigcapture_' scenario],'S','Sc','Mlist');
xx = 1:length(S);
xc = 1:length(Sc);
% plot(xx,1-cumsum(S),xc,1-cumsum(Sc));
% plot(xx,1-cumsum(S),xc,1-cumsum(Sc));
plot(xx,1-cumsum(S),xx,1-xx*S(1),'--');
xlabel('n_r');
ylabel('mean-square attenuation');
legend({'attenuation','approximation'});
clipboard('copy',['sigma_capture_' scenario]);
axis([1 60 0 1])
return;
