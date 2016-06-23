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

load C:\Dropbox\Osci\DATA\WeddKyst\hydrography/CTD_hist.mat
%%
clear ctdd

addpath('C:\Dropbox\Osci\FISP\#DATA\CTD')
for ipi = 1:numel(stns)
    figure(1);clf
    
    for i =1:numel(CTD)
        plot(CTD(i).Salt,CTD(i).PTemp,'.','color',0.8*[1 1 1])
        hold on
    end
    stn = stns{ipi};
    use(stn)
    switch stn
        case 'FSW1'
            cn = 1:5;
            col = [0*[1 1 1];0.3*[1 1 1]];
        case 'FSE2'
            cn = 1:6;
            col = [0.3*[1 1 1];0*[1 1 1]];
    end
    
    fl = 1;
    % plot(vfilt(mc.s,fl)',vfilt(mc.th,fl)','.','color','b')%;,cmap(n,:));
    x = vfilt(mc.s(:,cn),fl)';
    y = vfilt(mc.th(:,cn),fl)';
    z = mod(vfilt(yearfrac(mc.num(:,cn)),fl)',1)*360;
    scatter(x(:),y(:),40,z(:))%;,cmap(n,:));
    colormap(flipud(paruly(64)))
    hold on
   
    fl = 48;
    % plot(vfilt(mc.s,fl)',vfilt(mc.th,fl)','.','color','b')%;,cmap(n,:));
    x = vfilt(mc.s(:,cn),fl)';
    y = vfilt(mc.th(:,cn),fl)';
    z = mod(vfilt(yearfrac(mc.num(:,cn)),fl)',1)*360;
    %scatter(x(:),y(:),10,z(:),'filled')%;,cmap(n,:));
    plot(x',y','r','linewidth',1)
    
    cwin = [0 200];
    hold on
    
    odir = pwd;
    ctdd{1} = LoadFISSCTDData('FSW1');
    ctdd{2} = LoadFISSCTDData('FSE2');
    cd(odir);
    
    cs{1} = 5;
    cs{2} = 1:6;
    for i_ =1:numel(ctdd)
        ctd = ctdd{i_};
        clear s p th fla
        
    for i = cs{i_}
         dum = sw_salt(10*ctd(i).Cond/sw_c3515,ctd(i).Temp,ctd(i).Pres);
         dum(isnan(dum))=-inf;
         s{i} = dum;
         dum = ctd(i).Pres;
         dum(isnan(dum))=-inf;
         p{i} = dum;
         dum = sw_ptmp(s{i},ctd(i).Temp,p{i},0);
         dum(isnan(dum))=-inf;
         th{i} = dum;
         dum = ctd(i).Flag;
         dum(isnan(dum))=-inf;
         fla{i} = dum;
    end
    cell2col(s,p,th,fla);
    col2mat(s,p,th,fla);
    plot(s,th,'color',col(i_,:));
    end
    
    
    sref = 34.2:0.05:34.9;
    %D{n} = num2str(round(median(P)));
    plot(sref,...
        sw_ptmp(sref,...
        fp_t(sref,ones(size(sref))*vmedian(mc.p(:,6),1)),...
        ones(size(sref))*vmedian(mc.p(:,6),1),0),...
        '-','linewidth',2,'color',0.6*[0 1 0])
    
%      plot(sref,...
%           fp_t(sref,ones(size(sref))*vmedian(mc.p(:,6),1)),...
%                 '-','linewidth',2,'color',0.6*[1 1 1])
    plot(sref,thgade(sref,34.49,-2.45,'temp'),'-b')
    plot(sref,thgade(sref,34.51,-2.5,'temp'),'-r')
    box on
    %xlabel('Date')
    %ylabel(tit)
    
    
end

colorbar
caxis(cwin)
xlabel('Salinity')
ylabel('Potential temperature')
%legend(pl,D,'location','northwest')
title(stn)

for i=1:2
    if i==1
        xlim([34.2, 34.9]);
        ylim([-2.5146,-1.5]);
    elseif i==2
        xlim([34.4671, 34.6328]);
        ylim([-2.5146,-2.2248]);
    end
denscont(0);

thscr2png(['af_TS_' stn '_' datestr(now,'yyyy_mm') '_' num2str(i)],'150',outpath)
end
%%
