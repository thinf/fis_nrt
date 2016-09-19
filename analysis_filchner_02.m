% %% stack Microcat data
% clear num p t c s th
% 
% for n = 1:numel(data(1).Microcats)
%     num{n} = [];
%     p{n} = [];
%     t{n} = [];
%     c{n} = [];
%     
%     for i = 1:numel(data)
%         hold on
%         num_ = cell2mat(data(i).Microcats(n).Timestamp);
%         
%         p_ = data(i).Microcats(n).Pressure;
%         t_ = data(i).Microcats(n).Temperature;
%         c_ = data(i).Microcats(n).Conductivity;
%         
%         %abs(datenum(2016,1,1)-num)<500
%         ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<500 & p_ > 0 & p_ < 10000);
%         
%         num{n}=[num{n}; num_(ii)];
%         p{n}=[p{n}; p_(ii)];
%         t{n}=[t{n}; t_(ii)];
%         c{n}=[c{n}; c_(ii)];
%         
%         if NumberOfDailyMessages==5
%             num{n}(num{n}>datenum(2016,1,21))=datenum(2016,1,21);
%         end
%     end
%     
%     
%     %make dum p t th c s
%     %mc{n} = dum;
% end
% cell2col(num,p,t,c)
% col2mat(num,p,t,c)
% s = sw_salt(10*c/sw_c3515,t,p);
% th = sw_ptmp(s,t,p,0);
% 
% make mc num p t th c s
% 
% %% stack Aquadopp data
% clear num p t cv w head pitch roll amp1 amp2 amp3
% 
% for n = 1:numel(data(1).Aquadopps)
%     num{n} = [];
%     p{n} = [];
%     t{n} = [];
%     cv{n} = [];
%     w{n} = [];
%     head{n} = [];
%     pitch{n} = [];
%     roll{n} = [];
%     amp1{n} = [];
%     amp2{n} = [];
%     amp3{n} = [];
%     
%     for i = 1:numel(data)
%         hold on
%         %             for j = 1:numel(data(i).Aquadopps(n).Timestamp)
%         %                 if isempty(data(i).Aquadopps(n).Timestamp{j})
%         %                     data(i).Aquadopps(n).Timestamp{j} = nan;
%         %                 end
%         %             end
%         if isempty( data(i).Aquadopps(n).Timestamp{1})
%             data(i).Aquadopps(n).Timestamp{1} =nan;
%         end
%         num_ = data(i).Aquadopps(n).Timestamp{1}+[0:2:22]'/24; %cell2mat(data(i).Aquadopps(n).Timestamp);
%         
%         p_ = data(i).Aquadopps(n).P;
%         t_ = data(i).Aquadopps(n).T;
%         cv_ = data(i).Aquadopps(n).U + 1i*data(i).Aquadopps(n).V;
%         w_ = data(i).Aquadopps(n).W;
%         head_ = data(i).Aquadopps(n).Head;
%         pitch_ = data(i).Aquadopps(n).Pitch;
%         roll_ = data(i).Aquadopps(n).Roll;
%         amp1_ = data(i).Aquadopps(n).Amp1;
%         amp2_ = data(i).Aquadopps(n).Amp2;
%         amp3_ = data(i).Aquadopps(n).Amp3;
%         
%         %abs(datenum(2016,1,1)-num)<500
%         ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<500 & abs(cv_)<1);
%         
%         num{n}=[num{n}; num_(ii)];
%         p{n}=[p{n}; p_(ii)];
%         t{n}=[t{n}; t_(ii)];
%         cv{n}=[cv{n}; cv_(ii)];
%         w{n}=[w{n}; w_(ii)];
%         head{n}=[head{n}; head_(ii)];
%         pitch{n}=[pitch{n}; pitch_(ii)];
%         roll{n}=[roll{n}; roll_(ii)];
%         amp1{n}=[amp1{n}; amp1_(ii)];
%         amp2{n}=[amp2{n}; amp2_(ii)];
%         amp3{n}=[amp3{n}; amp3_(ii)];
%         
%         if NumberOfDailyMessages==5
%             num{n}(num{n}>datenum(2016,1,21))=datenum(2016,1,21);
%         end
%     end
%     
%     
%     %make dum p t th c s
%     %mc{n} = dum;
% end
% cell2col(num,p,t,cv,w,head,pitch,roll,amp1,amp2,amp3)
% col2mat(num, p,t,cv,w,head,pitch,roll,amp1,amp2,amp3)
% 
% make ad num p t cv w head pitch roll amp1 amp2 amp3

 stns{1} = 'FSW1';
 %stns{2} = 'FSE2';
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
for ipi = 1:numel(stns)
    stn = stns{ipi};
%% simple plot
%Not the prettiest way to compute variance ellipses
use(stn)
use ad
for i=1:size(cv,2)
    ii = find(isfinite(cv(:,i)));
    covxx(i,1)=vmean((real(cv(ii,i))-vmean(real(cv(ii,i)),1)).^2,1);
    covyy(i,1)=vmean((imag(cv(ii,i))-vmean(imag(cv(ii,i)),1)).^2,1);
    covxy(i,1)=vmean((real(cv(ii,i))-vmean(real(cv(ii,i)),1)).*(imag(cv(ii,i))-vmean(imag(cv(ii,i)),1)),1);
    cvm(i,1) = vmean(cv(ii,i),1);
end
depths = 1:4;
%Diagonalize the covariance matrix
[d1,d2,th]=specdiag(covxx,covyy,covxy);
%Convert to ellipse parameters, as in notes
a=sqrt(d1);b=real(sqrt(d2));  %bn can sometimes have a very small spurious negative part
[kappa,lambda]=ab2kl(a,b);  %Convert axes length a and b to RMS amplitude and 'linearity'

if strcmp(stn,'FSW1')
    kappa(1) = nan;
    th(1) = nan;
    lambda(1) = nan;
    cvm(1)=nan;
    j0 = 2;
else
    j0=1;
end
xwin = [-0.4 0.4];
ywin = [-0.4 0.4];
figure(500);clf
for i=j0:size(cv,2)
    subplot(2,2,i)
    plot(cv(:,i),'.','color',0.7*[1 1 1])
    hold on
    ellipseplot(kappa(i),lambda(i),th(i),cvm(i))
    plot([0 cvm(i)],'linewidth',2)
    xlim(xwin)
    ylim(ywin)
    hlines(0,'--k')
    vlines(0,'--k')
end
thscr2png(['velscat_' stn],'150',outpath)

xwin = [-0.2 0.2];
ywin = xwin;
figure(29);clf
cmap = jet(4);

hold on
for j = j0:4
    eh(j) = ellipseplot(kappa(j),lambda(j),th(j));
    set(eh(j),'color',cmap(j,:),'linewidth',2) % cahange ellipse color
    ph(j) = plot([0 real(cvm(j))], [0 imag(cvm(j))],'color',cmap(j,:),'linewidth',2); % add mean vector
end
title(stn)
xlim(xwin)
ylim(ywin)
hlines(0,'--k')
vlines(0,'--k')

lgd = ['1';'2';'3';'4'];
legend(ph(j0:end),lgd(j0:end),'location','northwest')
thscr2png(['var_ell_' stn],'150',outpath)
%Variance ellipses become larger and more circular near surface
%Also more angled north-south towards surface
%Orientation agrees with excess of downstream energy

figure(28);clf
subplot(1,3,1),plot(kappa,depths),title('Ellipse amplitude'),flipy
subplot(1,3,2),plot(lambda,depths),title('Ellipse linearity'),flipy
subplot(1,3,3),plot(th*(360/2/pi),depths); hold on; plot(angle(cvm)*(360/2/pi),depths,'r');
title('Ellipse orientation'),flipy
packcols(1,3)
thscr2png(['ellshape_' stn],'150',outpath)
%%

figure(400+ipi);clf
plot(ad.num,abs(ad.cv));yoffset(1);
title(stn)
ylabel('speed (m/s)')
datetick('x')
%xlim([datenum(2016,1,5), max(X)])
% xlim([min(ad.num), max(ad.num)])
%     grid on
box on
xlabel('Date')
thscr2png(['speed_ser_' stn],'150',outpath)

% simple plot
figure(200);
clf;for i = 1:4;hold on;ii = find(isfinite(ad.num(:,i)));scatter(ad.num(ii,i),abs(ad.cv(ii,i)),5,ad.t(ii,i)-median(ad.t(ii,i)),'filled');caxis([-0.02 0.02]);end;yoffset(1)
end

%% backscatter
figure(2020);clf
lgd = ['1';'2';'3';'4'];
for ipi = 1:numel(stns)
    stn = stns{ipi};

use(stn)
amp = ad.amp1+ad.amp2+ad.amp3;
subplot(1,2,ipi)
fl = 10;
plot(ad.num,vfilt(amp,fl,'median'))
ylim([100 250])
legend(lgd,'location','best')
title(['AD Backscatter, ' stn])
xlim([datenum(2015,12,24) now])
datetick('x')
end
thscr2png(['Backscatter'],'150',outpath)