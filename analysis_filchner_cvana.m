clear
close all
cd('C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\fis_nrt')

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

% coef(c)=ut_solv(mtime,real(cv),imag(cv),lat,constit(c),...
%        'OLS','white','NoDiagn','RunTimeDisp','nnn');
for ipi = 1:numel(stns)
    
    
    stn = stns{ipi};
    use(stn)
    switch stn
        case 'FSW1'
            cn = 2:4;
            ccn=2:6;
            
        case 'FSE2'
            cn = 1:4;
            ccn = 1:6;
            
    end
    %%
    figure(1);clf
    xwin=[datenum(2016,1,15) now];
    
    
    subplot(5,1,1:4)
    ywin = [-2.5 -2.2];
    plot(mc.num(:,ccn),mc.th(:,ccn))
    
    grid on
    xlim(xwin)
    ylim(ywin)
    datetick('x','keeplimits')
    set(gca,'XAxisLocation','top')
    
    subplot(5,1,5)
    ywin=[0 0.5];
    plot(ad.num,abs(ad.cv),'.','color',0.7*[1 1 1])
    hold on
    plot(ad.num,vfilt(abs(ad.cv),24))
    
    grid on
    xlim(xwin)
    ylim(ywin)
    datetick('x','keeplimits')
    
    thscr2png(['cvana_tser_' stn '_' datestr(now,'yyyy_mm')],'150',outpath)
    %%
    figure(2);clf
    
    subplot(5,1,1:4)
    xwin=[datenum(2016,1,15) now];
    ywin = [34.51 34.64];
    plot(mc.num(:,ccn),mc.s(:,ccn))
    
    grid on
    xlim(xwin)
    %ylim(ywin)
    datetick('x','keeplimits')
    set(gca,'XAxisLocation','top')
    
    subplot(5,1,5)
    ywin=[0 0.5];
    plot(ad.num,abs(ad.cv),'.','color',0.7*[1 1 1])
    hold on
    plot(ad.num,vfilt(abs(ad.cv),24))
    
    grid on
    xlim(xwin)
    ylim(ywin)
    datetick('x','keeplimits')
    thscr2png(['cvana_sser_' stn '_' datestr(now,'yyyy_mm')],'150',outpath)
end

%%
if 0
%%
    clear x y y_ x_
    x = nans(size(ad.cv(:,cn)));
   y_ = (vfilt((ad.cv(:,cn)),fl));
   y = nans(size(x));
    for j = 1:numel(cn)
        [x_, is] = sort(ad.num(:,cn(j)));
        y(:,j) = y_(is,j);
        [dum ia ic] = unique(x_);
        x(:,j)=nan;
        x(1:numel(dum),j)=dum;
        dum = y(ia,j);
        y(:,j)=nan;
        y(1:numel(dum),j)=dum;
    end
    
    
    %%
    % coef(c)=ut_solv(mtime,real(cv),imag(cv),lat,constit(c),...
%        'OLS','white','NoDiagn','RunTimeDisp','nnn');
constit={'M2','S2','N2','K1'};
    lat = -81.07583;
    %for i=1:size(x,2)
        %ii= find(isfinite(y(:,i)) & isfinite(x(:,i)));
    coef= ut_solv(x,real(y),imag(y),...
        lat,constit,'OLS','white','LinCI');%;,'white','NoDiagn','RunTimeDisp','nnn');
    %end
    
    [ uf, vf ] = ut_reconstr (x, coef );
    cvf = uf+1i*vf;
 %%
 figure(1);clf
ccn= [6 5 2 1];
xbin = [-0.4:0.04:0.4];
ybin = xbin;

 for i = 1:numel(cn)
     xmc = mc.num(:,(ccn(cn(i))));
     ymc = mc.p(:,(ccn(cn(i))));
     ii = find(isfinite(xmc) & isfinite(ymc));
     
     %ivar = interp1(xmc(ii),ymc(ii),x(:,i));
     ivar = interp1(xmc(ii),ymc(ii),x(:,i));
     ivar(1:30)=nan;
     ivar = ivar -vmedian(ivar,1);
     cwin = [-7*vstd(ivar,1) +7*vstd(ivar,1)];
     
     
     
 [mz,xmid,ymid,numz,stdz]= twodstats(real(y(:,i)),imag(y(:,i)),...
     ivar,...
     xbin,ybin);
 subplot(numel(cn),4,4*i-3)
 pcolor(xmid,ymid,mz); shading flat
 axis square
 colorbar
 caxis(cwin)
 hlines(0,'--k')
 vlines(0,'--k')
 text(-0.3,0.3,'p')
 end
 
  for i = 1:numel(cn)
     xmc = mc.num(:,(ccn(cn(i))));
     ymc = mc.th(:,(ccn(cn(i))));
     ii = find(isfinite(xmc) & isfinite(ymc));
     
     %ivar = interp1(xmc(ii),ymc(ii),x(:,i));
     ivar = interp1(xmc(ii),ymc(ii),x(:,i));
     ivar(1:30)=nan;
     cwin = [vmedian(ivar,1)-2*vstd(ivar,1) vmedian(ivar,1)+2*vstd(ivar,1)];
     
     
     
 [mz,xmid,ymid,numz,stdz]= twodstats(real(y(:,i)),imag(y(:,i)),...
     ivar,...
     xbin,ybin);
 subplot(numel(cn),4,4*i-2)
 pcolor(xmid,ymid,mz); shading flat
 axis square
 colorbar
 caxis(cwin)
 hlines(0,'--k')
 vlines(0,'--k')
 text(-0.3,0.3,'t')
 end

 %%
    figure(1);clf
    subplot(1,2,1)
    %z = mod(vfilt(yearfrac(mc.num(:,cn)),fl)',1)*360;
    plot(y,'.')%;,cmap(n,:));
    yoffset 0.4
    %plot([[0 0 0 0]; coef.umean],[[0 0 0 0]; coef.vmean])
    %colormap(flipud(paruly(64)))
    hold on
   
end
% 
% colorbar
% caxis(cwin)
% xlabel('Salinity')
% ylabel('Potential temperature')
% %legend(pl,D,'location','northwest')
% title(stn)
% 
% for i=1:2
%     if i==1
%         xlim([34.2, 34.9]);
%         ylim([-2.5146,-1.5]);
%     elseif i==2
%         xlim([34.4671, 34.6328]);
%         ylim([-2.5146,-2.2248]);
%     end
% denscont(0);
% 
% thscr2png(['af_TS_' stn '_' datestr(now,'yyyy_mm') '_' num2str(i)],'150',outpath)
% end
%%
