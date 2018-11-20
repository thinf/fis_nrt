 %% stack Microcat data
%  addpath('./nrt_toolbox')
    clear num p t c s th mc
    
    for n = 1:numel(data(1).Microcats)
        num{n} = [];
        p{n} = [];
        t{n} = [];
        c{n} = [];
        
        for i = 1:numel(data)
           
            num_ = cell2mat(data(i).Microcats(n).Timestamp);
            
            p_ = data(i).Microcats(n).Pressure;
            t_ = data(i).Microcats(n).Temperature;
            c_ = data(i).Microcats(n).Conductivity;
            
            %abs(datenum(2016,1,1)-num)<500
            ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<360*5.5 & p_ > 0 & p_ < 10000);
            
            num{n}=[num{n}; num_(ii)];
            p{n}=[p{n}; p_(ii)];
            t{n}=[t{n}; t_(ii)];
            c{n}=[c{n}; c_(ii)];
            
            if NumberOfDailyMessages==5
                num{n}(num{n}>datenum(2016,1,21))=datenum(2016,1,21);
            end
        end
        
        
        %make dum p t th c s
        %mc{n} = dum;
    end
    %%
%     cell2col(num,p,t,c)
%     col2mat(num,p,t,c)
%     s = sw_salt(10*c/sw_c3515,t,p);
%     th = sw_ptmp(s,t,p,0);
%     
%     make mc num p t th c s

mc.num = num;
mc.t = t;
mc.c = c;
mc.p = p;
    %% stack Aquadopp data
    clear num p t cv w head pitch roll amp1 amp2 amp3 ad
    
    for n = 1:numel(data(1).Aquadopps)
        num{n} = [];
        p{n} = [];
        t{n} = [];
        cv{n} = [];
        w{n} = [];
        head{n} = [];
        pitch{n} = [];
        roll{n} = [];
        amp1{n} = [];
        amp2{n} = [];
        amp3{n} = [];
        
        for i = 1:numel(data)
            hold on
            %             for j = 1:numel(data(i).Aquadopps(n).Timestamp)
            %                 if isempty(data(i).Aquadopps(n).Timestamp{j})
            %                     data(i).Aquadopps(n).Timestamp{j} = nan;
            %                 end
            %             end
            if isempty( data(i).Aquadopps(n).Timestamp{1})
                data(i).Aquadopps(n).Timestamp{1} =nan;
            end
            num_ = data(i).Aquadopps(n).Timestamp{1}+[0:2:22]'/24; %cell2mat(data(i).Aquadopps(n).Timestamp);
            
            p_ = data(i).Aquadopps(n).P;
            t_ = data(i).Aquadopps(n).T;
            cv_ = data(i).Aquadopps(n).U + 1i*data(i).Aquadopps(n).V;
            w_ = data(i).Aquadopps(n).W;
            head_ = data(i).Aquadopps(n).Head;
            pitch_ = data(i).Aquadopps(n).Pitch;
            roll_ = data(i).Aquadopps(n).Roll;
            amp1_ = data(i).Aquadopps(n).Amp1;
            amp2_ = data(i).Aquadopps(n).Amp2;
            amp3_ = data(i).Aquadopps(n).Amp3;
            
            %abs(datenum(2016,1,1)-num)<500
            ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<360*5.5 & abs(cv_)<1);
            
            num{n}=[num{n}; num_(ii)];
            p{n}=[p{n}; p_(ii)];
            t{n}=[t{n}; t_(ii)];
            cv{n}=[cv{n}; cv_(ii)];
            w{n}=[w{n}; w_(ii)];
            head{n}=[head{n}; head_(ii)];
            pitch{n}=[pitch{n}; pitch_(ii)];
            roll{n}=[roll{n}; roll_(ii)];
            amp1{n}=[amp1{n}; amp1_(ii)];
            amp2{n}=[amp2{n}; amp2_(ii)];
            amp3{n}=[amp3{n}; amp3_(ii)];
            
            if NumberOfDailyMessages==5
                num{n}(num{n}>datenum(2016,1,21))=datenum(2016,1,21);
            end
        end
        
        
        %make dum p t th c s
        %mc{n} = dum;
    end
    %%
%     cell2col(num,p,t,cv,w,head,pitch,roll,amp1,amp2,amp3)
%     col2mat(num, p,t,cv,w,head,pitch,roll,amp1,amp2,amp3)
%     
%     make ad num p t cv w head pitch roll amp1 amp2 amp3
%     timestamp = now;
%     make dum mc ad timestamp files path
%     eval([stn '= dum;'])
%     save([outpath stn '_' datestr(now,'yyyy_mm')],stn)
ad.num = num;
ad.cv = cv;

    %######## END OF DATASAVING SECTION, make standard plots ###############


%%
figure(1); set(gcf,'position',[100 10 1200 800]); clf

%
if strcmp(stn,'fsw1') || strcmp(stn,'fse2')
 adp =[6 5 2 1];
 numu = [datenum(2016,01,06):1:round(now)-1];

elseif strcmp(stn,'fne1') || strcmp(stn,'fne2')
 adp =[5 2 1];
 numu = [datenum(2016,12,01):1:round(now)-1];

end

clear tu pu

if strcmp(stn(1:4),'fsw1')
    ni = 1:numel(mc.num)-1;
else
    ni = 1:numel(mc.num);
end

for i = ni
    ii = i + numel(mc.num)-numel(ni);
    [c ia ic] = unique(mc.num{ii});
    %tu(:,i) = interp1(mc.num{ii},vfilt_loc(mc.t{ii},12),numu);
    %pu(i) = median(mc.p{ii});
    tu(:,i) = interp1(c,vfilt_loc(mc.t{ii}(ia),12),numu);
    pu(i) = median(mc.p{ii}(ia));
end
sp1 = subplot(4,1,1:2);
contourf(numu,pu,tu',[-2.5:0.025:-1.9],'color','none'); shading flat;
hold on
plot(numu(1),pu,'k>','linewidth',2,'markersize',10)
plot(numu(end),pu,'k<','linewidth',2,'markersize',10)
set(gca,'ydir','rev','xaxislocation','top','fontsize',16)
ylabel('Pressure [dBar]')
%title('Time')
datetick('x','keeplimits')
 cb = colorbar;
 set(cb,'fontsize',16)
 ylabel(cb,'Temperature [^oC]')
 caxis([-2.5 -1.9])
 colormap(flipud(hsv))



sp2=subplot(4,1,3);
sp1pos =  get(sp1,'position');
sp2pos =  get(sp2,'position');
sp2pos(3) = sp1pos(3)*0.925;


 clear adpu
 for i = 1:numel(ad.num)
    adpu{i} = [num2str(round(median(mc.p{adp(i)}))) ' dBar'];
 end


 if strcmp(stn,'fsw1') || strcmp(stn,'fse2')
  p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
        ad.num{2},abs(ad.cv{2})*100,'.',...
        ad.num{3},abs(ad.cv{3})*100,'.',...
        ad.num{4},abs(ad.cv{4})*100,'.');
    
    hold on
    ax=gca;ax.ColorOrderIndex=1;
    fl = 60;
    p2 = plot(ad.num{1},vfilt_loc(abs(ad.cv{1})*100,fl),...
        ad.num{2},vfilt_loc(abs(ad.cv{2})*100,fl),...
        ad.num{3},vfilt_loc(abs(ad.cv{3})*100,fl),...
        ad.num{4},vfilt_loc(abs(ad.cv{4})*100,fl),'linewidth',3);

elseif strcmp(stn,'fne1') || strcmp(stn,'fne2')
 p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
        ad.num{2},abs(ad.cv{2})*100,'.',...
        ad.num{3},abs(ad.cv{3})*100,'.');
    
    hold on
    ax=gca;ax.ColorOrderIndex=1;
    fl = 60;
    p2 = plot(ad.num{1},vfilt_loc(abs(ad.cv{1})*100,fl),...
        ad.num{2},vfilt_loc(abs(ad.cv{2})*100,fl),...
        ad.num{3},vfilt_loc(abs(ad.cv{3})*100,fl),'linewidth',3);
end
% 
%     p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
%         ad.num{2},abs(ad.cv{2})*100,'.',...
%         ad.num{3},abs(ad.cv{3})*100,'.',...
%         ad.num{4},abs(ad.cv{4})*100,'.');
%     
%     hold on
%     fl = 60;
%     p2 = plot(ad.num{1},vfilt_loc(abs(ad.cv{1})*100,fl),...
%         ad.num{2},vfilt_loc(abs(ad.cv{2})*100,fl),...
%         ad.num{3},vfilt_loc(abs(ad.cv{3})*100,fl),...
%         ad.num{4},vfilt_loc(abs(ad.cv{4})*100,fl),'linewidth',3);
    set(gca,'ylim',[0 50])
    ylabel('Speed [cm/s]','fontsize',16)
    set(sp2,'xlim',get(sp1,'xlim'),'fontsize',16)
    grid on
    datetick('x','keeplimits')
    set(sp2,'xticklabel','')
    lgd = legend(p2,adpu,'location','eastoutside');
set(sp2,'position',sp2pos)

%
sp3=subplot(4,1,4);
sp1pos =  get(sp1,'position');
sp3pos =  get(sp3,'position');
sp3pos(3) = sp1pos(3)*0.925;


 clear adpu
 for i = 1:numel(ad.num)
    adpu{i} = [num2str(round(median(mc.p{adp(i)}))) ' dBar'];
 end

fl = 12*30;
 if strcmp(stn,'fsw1') || strcmp(stn,'fse2')
%   p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
%         ad.num{2},abs(ad.cv{2})*100,'.',...
%         ad.num{3},abs(ad.cv{3})*100,'.',...
%         ad.num{4},abs(ad.cv{4})*100,'.');
    p1 = plot(ad.num{1},vfilt_loc(real(ad.cv{1})*100,fl),...
        ad.num{2},vfilt_loc(real(ad.cv{2})*100,fl),...
        ad.num{3},vfilt_loc(real(ad.cv{3})*100,fl),...
        ad.num{4},vfilt_loc(real(ad.cv{4})*100,fl),'linewidth',2);
    hold on
    ax=gca;ax.ColorOrderIndex=1;
    
    p2 = plot(ad.num{1},vfilt_loc(imag(ad.cv{1})*100,fl),'--',...
        ad.num{2},vfilt_loc(imag(ad.cv{2})*100,fl),'--',...
        ad.num{3},vfilt_loc(imag(ad.cv{3})*100,fl),'--',...
        ad.num{4},vfilt_loc(imag(ad.cv{4})*100,fl),'--','linewidth',2);

elseif strcmp(stn,'fne1') || strcmp(stn,'fne2')
%  p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
%         ad.num{2},abs(ad.cv{2})*100,'.',...
%         ad.num{3},abs(ad.cv{3})*100,'.');
%     
 p1 = plot(ad.num{1},vfilt_loc(real(ad.cv{1})*100,fl),'--',...
        ad.num{2},vfilt_loc(real(ad.cv{2})*100,fl),'--',...
        ad.num{3},vfilt_loc(real(ad.cv{3})*100,fl),'--','linewidth',3);   
hold on
    p2 = plot(ad.num{1},vfilt_loc(imag(ad.cv{1})*100,fl),...
        ad.num{2},vfilt_loc(imag(ad.cv{2})*100,fl),...
        ad.num{3},vfilt_loc(imag(ad.cv{3})*100,fl),'linewidth',3);
end
% 
%     p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
%         ad.num{2},abs(ad.cv{2})*100,'.',...
%         ad.num{3},abs(ad.cv{3})*100,'.',...
%         ad.num{4},abs(ad.cv{4})*100,'.');
%     
%     hold on
%     fl = 60;
%     p2 = plot(ad.num{1},vfilt_loc(abs(ad.cv{1})*100,fl),...
%         ad.num{2},vfilt_loc(abs(ad.cv{2})*100,fl),...
%         ad.num{3},vfilt_loc(abs(ad.cv{3})*100,fl),...
%         ad.num{4},vfilt_loc(abs(ad.cv{4})*100,fl),'linewidth',3);
    set(gca,'ylim',[-8 8])
    plot(get(sp1,'xlim'),[0 0],'k')
    ylabel('Speed [cm/s]','fontsize',16)
    set(sp3,'xlim',get(sp1,'xlim'),'fontsize',16)
    grid on
    datetick('x','keeplimits')
    lgd = legend([p1(1),p2(1)],'eastward','northward','location','eastoutside');
set(sp3,'position',sp3pos)
%error('stop')
%%
oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)

print('-dpng', [savepath '/' stn '_timeseries'], '-r300');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)
