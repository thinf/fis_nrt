%%
clear all
close all
fgpath = './working_plots/';
fgpf = 'af_tide_spectra_';
 stns{1} = 'FSW1';
 stns{2} = 'FSE2';
    % specify folder where temporary files will be created. Also the files
    % "blankmsg1.sbd" to "blankmsg7.sbd" MUST be stored here.
outpath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\processed/';

%%
for ipi = 1:numel(stns)
    stn = stns{ipi};
    
    try
        load([outpath stn '_' datestr(now,'yyyy_mm')])
    catch
        disp('no recent dataset found, use data from march')
        load([outpath stn '_' datestr(now,'yyyy_03')])
    end
end
%%
use FSE2
figure(1);clf
plot(mc.num,mc.t); legend('1','2','3','4','5','6')
ylim([-2.5 -2.15])
fgnme = 'temp_series';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])


figure(101);clf
plot(mc.p,mc.t,'.')
ylim([-2.5 -2.15])
xlim([720 1200])
fgnme = 'temp_vs_pressure';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

figure(102);clf
ad.num(ad.num-ad.num_logger<-0.9)=ad.num(ad.num-ad.num_logger<-0.9)+1;
plot(ad.num,abs(ad.cv))
fgnme = 'speed_series';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])


%%
% Tore,
% I've been staring at the data- this is some complex material.
% There're many interesting features- most of them you already pointed to.
% The currents seem quite barotropic, and some depth dependency in the
% barotropic tidal currents can be expected given the boundary layer and latitude effects.
% I doubt that there's much internal waves. Having said that, there seems
% to be +/- 20/50 m or so vertical displacement of isotherms from the CTD profiles
% in the stratified portion.
%
% Lateral gradient is strong. Spring-neap cycle is strong. Boundary layers are thick.
% It would be important to get hold of the cross-cavity bathymetry (bottom and ice).
% I suspect there can be some focussing/resonance of internal tide energy, if any.
% 
% I would start off with some spectral analysis.
% I would then attempt to calculate vertical isotherm displacements over 24-h or so
% moving average background dT/dz.

%%
% simple Lomb-Scargle Periodogram for non-uniformly spaced data
if 0
cmap=jet(6);
figure(1);clf
for i = 1:6;
    ii = find(isfinite(mc.num(:,i))&isfinite(mc.p(:,i)));
    [f,P,prob] = lomb(mc.num(ii,i)-mc.num(ii(1),i),mc.p(ii,i),12,2);
    
    % ii = find(isfinite(ad.num(:,i))&isfinite(abs(ad.cv(:,i))));
    %[f,P,prob] = lomb(ad.num(ii,i)-ad.num(ii(1),i),abs(ad.cv(ii,i)),12,2);
    
    subplot(2,1,1);hold on;plot(f,vfilt(P,24,'median'),'color',cmap(i,:));
    set(gca,'YScale','log','xscale','log')
    
    subplot(2,1,2);hold on;plot(f,vfilt(prob,24),'color',cmap(i,:));set(gca,'xscale','log')
    
end
% inertial waves:
T = 2*pi()/sw_f(83);

subplot(2,1,1)
vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
%vlines(bsxfun(@mtimes,[2 3 4]',tidefreq)*24,':k')
vlines([6],'2m')
%vlines([3 1.5],'m')

subplot(2,1,2)
vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
%vlines(bsxfun(@mtimes,[2 3 4]',tidefreq)*24,':k')
vlines([6],'2m')
%vlines([3 1.5],'m')

legend('2','3','4','5','6')
end
%% grid data
% FSW1 unified time vector:
%mc.num=mc.num(:,2:end);
%mc.t=mc.t(:,2:end);
%mc.p=mc.p(:,2:end);
 % numu = datenum(2016,1,23,18+[2:2:3.5*size(mc.num,1)],0,0);
  %pu = 750:2:1000;
  
use FSE2
% fix ad-time stamp
ad.num(ad.num-ad.num_logger<-0.9)=ad.num(ad.num-ad.num_logger<-0.9)+1;

% unified time stamp and pressure for microcats
numu = datenum(2016,1,8,0+[2:2:2.2*size(mc.num,1)],0,0);
pu = 700:2:1200;

% simplify time vector for microcats
 dv = datevec(mc.num(:));
 dv(:,5)=0;
 dv(:,6)=0;
 nums = datenum(dv);
 
 % initialize varaible and interpolate temperature data (vertical 1d)
 tu = nan(numel(numu),numel(pu));
 su = nan(numel(numu),numel(pu));
 pd = mc.p(:,1:end);
 td = mc.t(:,1:end);
 sd = mc.s(:,1:end);
 for i = 1:numel(numu)
     ii = find(numu(i)==nums);
     if numel(ii)>1
        % tu(i,:) = interp1(pd(ii),td(ii),pu);
        disp('NOTE: using cubic interpolation, watch out for temp inversion near ice base')
        tu(i,:) = interp1(pd(ii),td(ii),pu,'cubic',nan);
        su(i,:) = interp1(pd(ii),sd(ii),pu,'linear',nan);
     end
 end
% FIXME: Salinity / stratification (N2) profiles need Salinity offset correction

% mc.t(:,1)=nan;
%  ii = find(isfinite(mc.num) & isfinite(mc.p) & isfinite(mc.t));
%  tu = griddata(mc.num(ii),mc.p(ii)/100000,mc.t(ii),numu,pu'/100000)';
 
 % simplify time vector for ads
 dvc = datevec(ad.num(:));
 dvc(:,5)=0;
 dvc(:,6)=0;
 numcs = datenum(dvc);

% fix pressure coordinate for aquadopps
 pdc = repmat(vmedian(mc.p(:,[1 2 5 6]),1),[size(ad.cv,1),1]);
 
 % interpolate velocity data on uniform grid (2d interp with most weigthing
 % on vertical data)
 ii = find(isfinite(ad.num) & isfinite(pdc) & isfinite(ad.cv));
 cvu = griddata(ad.num(ii),pdc(ii)/100000,ad.cv(ii),numu,pu'/100000)';
 
 %%
 thfig(201);
 subplot(2,1,1)
  pcolor(numu,pu,vfilt(abs(cvu),48)');shading flat
  hlines(pdc(1,:))
  title('spring neap')
  set(gca,'ydir','rev')
  
  subplot(2,1,2)
  pcolor(numu,pu,bsxfun(@minus,vfilt(angle(cvu),48),vmean(vfilt(angle(cvu),48),2))');shading flat
  hold on
  contour(numu,pu,bsxfun(@minus,vfilt(angle(cvu),48),vmean(vfilt(angle(cvu),48),2))',[0 0],'k');shading flat
  hlines(pdc(1,:))
  title('vertical directional shear')
  set(gca,'ydir','rev')
  
fgnme = 'gridded_velocity';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])
%%
ipu = 130;
thfig(3)
ax(1) = subplot(3,1,1);
pcolor(numu,pu,tu')
shading flat
hold on
set(gca,'ydir','rev')
axis tight
plot(mc.num,mc.p,'k')

hlines(pu(ipu),'G')

ax(2) = subplot(3,1,2);
tuz = tu(:,ipu);
tmedu = bsxfun(@minus,tu,vmedian(tu,1));
%tfilt = vfilt(vfilt(tu(:,ipu),24,'median'),12);
%tfilt = vfilt(vfilt(tu(:,ipu),24,'median'),240);
tfilt = vfilt(vfilt(tu(:,ipu),48,'median'),240);
%plot(numu,tmedu(:,ipu))
plot(numu,tuz)
hold on
plot(numu,tfilt,'r')
%axis tight

ax(3) = subplot(3,1,3);
plot(numu,tuz-tfilt)
linkaxes(ax,'x')

fgnme = 'gridded_temperature';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])


%% lomb spectrum
if 0
ii = find(isfinite(numu')&isfinite(tuz)&isfinite(tfilt));
[f,P,prob] = lomb((numu(ii)-numu(1))',tuz(ii),12,2);

figure(4);clf
  
    subplot(2,1,1);hold on;plot(f,vfilt(P,24,'median'));
    set(gca,'YScale','log','xscale','log')
    
    subplot(2,1,2);hold on;plot(f,vfilt(prob,24));set(gca,'xscale','log')
  
% inertial waves:
T = 2*pi()/sw_f(83);

subplot(2,1,1)
vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
%vlines(bsxfun(@mtimes,[2 3 4]',tidefreq)*24,':k')
vlines([6],'2m')
%vlines([3 1.5],'m')

subplot(2,1,2)
vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
%vlines(bsxfun(@mtimes,[2 3 4]',tidefreq)*24,':k')
vlines([6],'2m')
%vlines([3 1.5],'m')

legend('2','3','4','5','6')
end

%% multitaper Temperature
tuzf = tuz;
tuzf(isnan(tuz))=tfilt(isnan(tuz));
ii = find(isfinite(numu')&isfinite(tuzf)&isfinite(tfilt));


P = 6;
psi=sleptap(numel(ii),P);

%

thfig(5)
% inertial waves:
T = 2*pi()/sw_f(83);


ax1(1) = subplot(2,1,1);
[f,s]=mspec(1/12,tuzf(ii)-tfilt(ii),psi);
plot(f,s)
 set(gca,'YScale','log','xscale','log')
 vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
vlines([6],'2m')
title('high-passed temperature spectrum')
xlabel('f (1/d)')
ylabel('multitaper spectral resonance')

ax1(2) = subplot(2,1,2);
[f,s]=mspec(1/12,tuzf(ii),psi);
plot(f,s)
 set(gca,'YScale','log','xscale','log')
 vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
vlines([6],'2m')
title('raw temperature spectrum')
xlabel('f (1/d)')
ylabel('multitaper spectral resonance')

linkaxes(ax1,'x')

fgnme = 'temperature_spectrum';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%% multitaper velocity
cvuz = cvu(:,ipu);
cvfilt = vfilt(vfilt((cvuz),12,'median'),48);
tuzf = (cvuz);
tuzf(isnan(cvuz))=(cvfilt(isnan(cvuz)));
ii = find(isfinite(numu')&isfinite(tuzf)&isfinite(cvfilt));

P = 4;
psi=sleptap(numel(ii),P);

%
figure(5);clf

ax1(1) = subplot(2,1,1);
[f,spp,snn,spn]=mspec(1/12,cvuz(ii),psi);
% [f,s]=mspec(1/12,abs(tuzf(ii))-abs(cvfilt(ii)),psi);
plot(f,spp,f,snn)
linestyle b r g
 set(gca,'YScale','log','xscale','log')
 vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
vlines([6],'2m')
title('complex rotary spectrum')
xlabel('f (1/d)')
ylabel('multitaper spectral resonance')

ax1(2) = subplot(2,1,2);
[f,s]=mspec(1/12,abs(tuzf(ii)),psi);
plot(f,s)
 set(gca,'YScale','log','xscale','log')
 vlines([1 2 3]./T*3600*24,'2G')
vlines(tidefreq*24,'--k')
vlines([6],'2m')
title('Absolute velocity spectrum')
xlabel('f (1/d)')
ylabel('multitaper spectral resonance')

linkaxes(ax1,'x')

fgnme = 'velocity_multitaper_spectrum';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%% wavelet transform

ga=3;be=1;
fs=morsespace(ga,be,{0.05,pi},2*pi/1200,8);
[psi,psif]=morsewave(numel(ii),1,ga,be,fs(end));
[wp,wn]=wavetrans(cvuz(ii),conj(cvuz(ii)),{1,3,1,fs,'bandpass'},'mirror');
[wp,wn]=wavetrans(cvuz(ii),conj(cvuz(ii)),{1,3,1,fs,'bandpass'},'mirror');
    

    thfig(501)
    h=wavespecplot(1:numel(ii),vfilt(cvuz(ii),12),fs,log(sqrt(abs(wp).^2+abs(wn).^2)));
   
    hlines(tidefreq('k1')*12), hlines(tidefreq('m2')*12),hlines(tidefreq('mf')*12)
    hlines(1./T*3600*12,'2G')
    caxis([-6 -1])
title('low passed currents and wavelet transform')

fgnme = 'velocity_wavelet_transform';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%% thermocline depth
tcls = [-2.25 -2.33 -2.4];
z_tcl = nans(numel(numu),numel(tcls));

for i = 1:numel(tcls)
    tcl = tcls(i);
tud = abs(tu-tcl); % look for thermocline at -2.3 deg
tud(tud>0.05)=nan;
[ii jj] = find(bsxfun(@minus,tud,min(tud,[],2))==0);

z_tcl(ii,i)=pu(jj);
end

ipu = 150;
thfig(6)

ax(1) = subplot(3,1,1);
pcolor(numu,pu,tu')
shading flat
hold on
set(gca,'ydir','rev')
axis tight
plot(mc.num,mc.p,'k')
plot(numu,z_tcl,'w','linewidth',1)
ylabel('pressure (dBar)')
title(['temperature profile and themocline depth at ' num2str(tcls) ' degC'])
colorbar('location','east')

%
ax(2) = subplot(3,1,2);
tclfilt = vfilt(vfilt(z_tcl,12,'median'),12);
%plot(numu,tmedu(:,ipu))
plot(numu,z_tcl)
hold on
plot(numu,tclfilt,'k')
%axis tight
title('raw and low pass filtered thermocline depth')
ylabel('pressure (dBar)')

ax(3) = subplot(3,1,3);
plot(numu,z_tcl-tclfilt)
linkaxes(ax,'x')
title('raw minus low pass filtered thermocline depth')
ylabel('pressure (dBar)')

fgnme = 'thermocline_depth';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%%
thfig(7)
% thermoclilne depth
axt(1) =subplot(3,1,1);
plot(numu,z_tcl(:,2))
hold on
plot(numu,vfilt(vfilt(z_tcl(:,2),12,'median'),48),'k','linewidth',2)
ylabel('thermocline depth')
axis tight
% ylim([57 177])
grid on
title('thermocline thickness correlates with speed; thermocline height lags speed and thickness by 3-4 days / is highest at decaying spring tide')

% thermoclilne thickness
axt(2) =subplot(3,1,2);
plot(numu,z_tcl(:,1)-z_tcl(:,3))
hold on
plot(numu,vfilt(vfilt(z_tcl(:,1)-z_tcl(:,3),12,'median'),48),'k','linewidth',2)
ylabel('thermocline thcikness')
axis tight
ylim([57 177])
grid on

axt(3) =subplot(3,1,3);
plot(numu,abs(cvuz));
hold on
plot(numu,vfilt(vfilt(abs(cvuz),12,'median'),48),'r','linewidth',2);
ylabel('velocity')
linkaxes(axt,'x')
ylim([0.005 0.37])
grid on

packrows(3,1)

fgnme = 'thermocline_depth_and_thickness_vs_velocity';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])


%% 2-D histograms
tcln=2;
thfig(8)

cvf = vfilt(vfilt((cvuz),12,'median'),48);
ztf = vfilt(vfilt(z_tcl,12,'median'),48);

[hm, xmid,ymid, h, hstd] = twodstats(real(cvuz-cvf),...
    imag(cvuz-cvf),-1*(z_tcl(:,tcln)-ztf(:,tcln)),...
    [-0.2:0.01:0.2],[-0.5:0.01:0.5]);

subplot(1,2,1)
pcolor(xmid,ymid,log(h));shading flat
vlines(0)
hlines(0)
legend('velocity histogram')

subplot(1,2,2)
pcolor(xmid,ymid,hm);shading flat
caxis([-50 50])
vlines(0)
hlines(0)
colorbar
legend('tcl depth anomaly')

fgnme = 'velocity_histogram_and_thermocline_anomly_distribution';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])


% 
% [hm, xmid,ymid, h, hstd] = twodstats(real(cvuz-cvf),...
%     imag(cvuz-cvf),z_tcl(:,1)-z_tcl(:,3),...
%     [-0.2:0.01:0.2],[-0.5:0.01:0.5]);
% 
% subplot(2,2,3)
% pcolor(xmid,ymid,log(h));shading flat
% vlines(0)
% hlines(0)
% subplot(2,2,4)
% pcolor(xmid,ymid,hm);shading flat
% %caxis([-50 50])
% vlines(0)
% hlines(0)
% colorbar
% legend('tcl depth')


%% low pass correaltions
thfig(9)

% a = vfilt(vfilt(abs(cvuz),12,'median'),48);
% b = vfilt(vfilt(z_tcl(:,1)-z_tcl(:,3),12,'median'),48);
% 
% i1 = max(find(isfinite(a),1,'first'),find(isfinite(b),1,'first'));
% i2 = min(find(~isfinite(a(i1:end)),1,'first'),find(~isfinite(b(i1:end)),1,'first'))+i1-2;
% 
%  disp('correlatio')
%  disp(corrcoef(a(i1:i2),b(i1:i2)))
% 

cvf = vfilt(vfilt(abs(cvuz),12,'median'),4*12);
tcldf = vfilt(vfilt(z_tcl(:,1)-z_tcl(:,3),12,'median'),4*12);
tclf = vfilt(vfilt(z_tcl(:,2),12,'median'),4*12);

i1 = max([find(isfinite(cvf),1,'first'),find(isfinite(tcldf),1,'first'),find(isfinite(tclf),1,'first')]);
i2 = min([find(~isfinite(cvf(i1:end)),1,'first'),find(~isfinite(tcldf(i1:end)),1,'first'),find(~isfinite(tclf(i1:end)),1,'first')])+i1-2;
ii = i1:i2;

subplot(2,3,1:3)
plot(numu(ii),thnorm(abs(cvf(ii))),numu(ii),thnorm((tcldf(ii))),numu(ii),thnorm((tclf(ii))))
legend('speed','thermocline thickness','thermocline depth','location','west')
title('low-passed correlations of temperature with currents')
xlabel('normalized amplitude')
ylabel('time (julian days)')
% left: b lags a

subplot(2,3,4)
maxlag = 28*12;
plot([-maxlag:maxlag]/12,xcorr(thnorm(abs(cvf(ii))),thnorm((tcldf(ii))),maxlag,'coeff'))
%vlines([maxlag, maxlag-3*12, maxlag+11*12])
vlines([0],'2k')
%vlines([-4*12, +11*12])
title('thickness is in phase with speed')
vlines(numel(ii))
xlabel('time lag (days)')
ylabel('xcorr(speed,thickness)')


subplot(2,3,5)
maxlag = 28*12;
plot([-maxlag:maxlag]/12,xcorr(thnorm(abs(cvf(ii))),thnorm((tclf(ii))),maxlag,'coeff'))
%vlines([maxlag, maxlag-3*12, maxlag+11*12])
vlines([0],'2k')
vlines([-3.5*12]/12)
title('thermocline height lags speed by 3.5 days')
vlines(numel(ii))
xlabel('time lag (days)')
ylabel('xcorr(speed,depth)')


subplot(2,3,6)
maxlag = 28*12;
plot([-maxlag:maxlag]/12,xcorr(thnorm((tcldf(ii))),thnorm((tclf(ii))),maxlag,'coeff'))
%vlines([maxlag, maxlag-3*12, maxlag+11*12])
vlines([0],'2k')
vlines([-3.5*12 9.5*12]/12)
title('thermocline height lags thermocline thickness by 3.5 days. Or leads by 9.5 days ?')
vlines(numel(ii))
xlabel('time lag (days)')
ylabel('xcorr(thickness,depth)')

fgnme = 'thermocline_temperature_correlations';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])


%% dpth distribution and associated circualtion modes
tcln = [2];
xbin = [-54:6:54];
xwin=[-60 60];

cvf = vfilt(vfilt((cvuz),12,'median'),48);
ztf = vfilt(vfilt(z_tcl,12,'median'),48);
    
tclhp = -1*(z_tcl(:,tcln)-ztf(:,tcln));
cvhp = (cvuz-cvf)*100;

% quantiles
x = tclhp(isfinite(tclhp));
% compute the median
medianx = median(x);

% STEP 1 - rank the data
y = sort(x);

% compute 25th percentile (first quartile)
Q(1) = median(y(find(y<median(y))));

% compute 50th percentile (second quartile)
Q(2) = median(y);

% compute 75th percentile (third quartile)
Q(3) = median(y(find(y>median(y))));

% compute Interquartile Range (IQR)
IQR = Q(3)-Q(1);

 
thfig(10)
ax1= subplot(2,1,1);


hist(tclhp,[-60:1:60])
%hist(-1*(z_tcl(:,tcln)-ztf(:,tcln)),xbin)
xlim(xwin)
 %vlines([Q(1) Q(3) Q(1)-IQR Q(3)+IQR],'-k')
 vlines([Q(1) Q(3)],'-k')

 
xlabel('thermocline depth anomaly (dBar, positive upward)')
ylabel('counts per bin')
title('low-pass filtered thermocline depth distribution, 2 dBar bins')

 %% variance ellipses
[hm, xmid,ymid, h, hstd] = twodstats(tclhp,ones(size(cvuz)),...
    cvhp,xbin,[0.5 1.5]);

covxx = twodstats(tclhp,ones(size(cvuz)),...
    (real(cvhp)-vmean(real(cvhp),1)).^2,...
    xbin,[0.5 1.5]);

covyy = twodstats(tclhp,ones(size(cvuz)),...
    (imag(cvhp)-vmean(imag(cvhp),1)).^2,...
    xbin,[0.5 1.5]);

covxy = twodstats(tclhp,ones(size(cvuz)),...
    (real(cvhp)-vmean(real(cvhp),1)).*(imag(cvhp)-vmean(imag(cvhp),1)),...
    xbin,[0.5 1.5]);


 % variance ellipses and mean vector

%Diagonalize the covariance matrix
[d1,d2,th]=specdiag(covxx,covyy,covxy);
%Convert to ellipse parameters, as in notes
a=sqrt(d1);b=real(sqrt(d2));  %bn can sometimes have a very small spurious negative part
[kappa,lambda]=ab2kl(a,b);  %Convert axes length a and b to RMS amplitude and 'linearity'

ax2 = subplot(2,1,2);

%thfig(12)
for i=1:numel(hm)
    
    plot(hm(i)+xmid(i),'o','color',0*[1 1 1])
   
    hold on
    %ellipseplot(kappa(i),lambda(i),gcf,hm(i)+xmid(i))
    el = ellipseplot(kappa(i),lambda(i),th(i),xmid(i));
    
    set(el,'color','k')
    plot([0 hm(i)]+xmid(i),'linewidth',2,'color','k')
    
    
end

ywin = [-20 20];
ylim(ywin)
  xlim(xwin)
    hlines(0,'--k')
    vlines(0,'--k')
     %vlines([Q(1) Q(3) Q(1)-IQR Q(3)+IQR],'-k')
     vlines([Q(1) Q(3)],'-k')
   
    set(ax2,'plotboxaspectratio',[diff(ywin),diff(xwin),1])
set(ax1,'plotboxaspectratio',get(ax2,'plotboxaspectratio'))
xlabel('easting speed (cm/s)')
ylabel('northing speed (cm/s)')
title('corresponding mean currents and variance ellipses, 6 dBar bins')

fgnme = 'thermocline_depth_distribution_and associated_currents';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%% hodograph
thfig(11)
 hodograph(hm,xmid)


%% compute modes
% mbins = sort([inf -inf Q(1) Q(3) Q(1)-IQR Q(3)+IQR]);
mbins = sort([inf -inf Q(1) Q(3)]);

[bnum xb bmid] = bindata(mbins,tclhp);

tuf = vfilt(vfilt(tu,12,'median'),48);
suf = vfilt(vfilt(su,12,'median'),48);
cvuf = vfilt(vfilt(cvu,12,'median'),48);

for i=1:numel(bmid)
    ii = find(bnum == i);
    tm(:,i) = vmean(tu(ii,:),1);
    sm(:,i) = vmean(su(ii,:),1);
    cvm(:,i) = vmean(cvu(ii,:),1);
    
    tmanom(:,i) = vmean(tu(ii,:)-tuf(ii,:),1);
    smanom(:,i) = vmean(su(ii,:)-suf(ii,:),1);
    cvmanom(:,i) = vmean(cvu(ii,:)-cvuf(ii,:),1);
end

thfig(12)
 ii = find(isfinite(tm));
hodograph(cvm(ii),tm(ii))


%% T-S
thfig(13)
plot(mc.s,mc.t,'.','color',0.7*[1 1 1])
hold on

pl = plot(sm,tm,'linewidth',2);
legend(pl,'low','normal','high','location','east')

    sref = 34.4:0.05:34.65;
    plot(sref,fp_t(sref,ones(size(sref))*vmedian(mc.p(:,end),1)),'-','linewidth',2,'color',0*[1 1 1])
    plot(sref,thgade(sref,34.49,-2.45,'temp'),'-b')
    plot(sref,thgade(sref,34.51,-2.5,'temp'),'-r')
    box on

xlim([34.4671, 34.6328]);
%ylim([-2.5146,-2.2248]); % potential temp
ylim([-2.5146,-2.12248]); % in-situ temp
denscont(0);
xlabel('Salinity')
ylabel('Potential temperature')

fgnme = 'TS_modes';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%% vertical modes
thfig(14)
subplot(1,3,1)
plot(tm,pu,'linewidth',2)
set(gca,'ydir','rev')
grid on
legend('low','normal','high','location','south')
ylabel('pressure (dBar)')
xlabel('Temperature (degC)')
title('Vertical temperature modes')

subplot(1,3,2)
plot(abs(cvm),pu,'linewidth',2)
set(gca,'ydir','rev')
grid on
ylabel('pressure (dBar)')
xlabel('speed (m/s)')
title('vertical velocity modes')

subplot(1,3,3)
plot(bsxfun(@minus,angle(cvm),vmean(angle(cvm),1))*180/pi(),pu,'linewidth',2)
set(gca,'ydir','rev')
vlines(0)
grid on
ylabel('pressure (dBar)')
xlabel('deviation from vertical mean direction (deg)')
title('vertical flow direction modes')

fgnme = 'vertical_modes';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%% vertical mode anomalies
thfig(15)

subplot(1,3,1)
plot(tmanom,pu,'linewidth',2)
set(gca,'ydir','rev')
grid on
legend('low','normal','high','location','south')
ylabel('pressure (dBar)')
xlabel('Temperature (degC)')
title('Vertical temperature mode anomalies')

subplot(1,3,2)
plot(abs(cvmanom),pu,'linewidth',2)
set(gca,'ydir','rev')
grid on
ylabel('pressure (dBar)')
xlabel('speed (m/s)')
title('vertical velocity mode anomalies')

subplot(1,3,3)
plot(bsxfun(@minus,angle(cvmanom),vmean(angle(cvmanom),1))*180/pi(),pu,'linewidth',2)
set(gca,'ydir','rev')
vlines(0)
grid on
ylabel('pressure (dBar)')
xlabel('deviation from vertical mean direction (deg)')
title('vertical flow direction mode anomalies')



fgnme = 'vertical_mode_anomalies';
thscr2png([fgpf fgnme],'150',fgpath)
saveas(gcf,[fgpath fgpf fgnme '.fig'])

%%
save('../FSE2_gridded','numu','pu','tu','cvu')