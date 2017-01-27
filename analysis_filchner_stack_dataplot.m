
%% plot Microcats
figure(1); clf

cmap = jet(6);

for d = 1:3
    subplot(3,1,d)
    
    for n = 1:6
        X = [];
        Y = [];
        for i = 1:numel(data)
            hold on
            x = cell2mat(data(i).Microcats(n).Timestamp);
            switch d
                case 1
                    
                    y = data(i).Microcats(n).Pressure;
                    ywin = [700 1200];
                    tit = 'Pressure (dBar)';
                    set(gca,'ydir','rev')
                case 2
                    %x = data(i).Microcats(n).Timestamp{:};
                    y = data(i).Microcats(n).Temperature;
                    ywin = [-2.5 -2.2];
                    tit = 'Temperature (degC)';
                case 3
                    %x = data(i).Microcats(n).Timestamp{:};
                    y = data(i).Microcats(n).Conductivity;
                    ywin = [2.68 2.75];
                    tit = 'Conductivity (mS/cm)';
            end
            abs(datenum(2016,1,1)-x)<200;
            ii = find(isfinite((x)) & abs(datenum(2016,1,1)-x)<300 & isfinite(y));
            X=[X; x(ii)];
            Y=[Y; y(ii)];
            
        end
        if NumberOfDailyMessages==5
            X(X>datenum(2016,1,21))=datenum(2016,1,21);
        end
        ii = find(abs(Y-median(Y))<4*std(Y));
        plot(X(ii),Y(ii),'-','color',cmap(n,:),'linewidth',2)
        %plot(X(ii),Y(ii),'-','color',cmap(n,:),'.')
        % tit = 'Pressure (dBar)';
        
        %     title(['AD ' num2str(n) ', color = time (days)'])
        %     xlabel(['U (cm/s)'])
        %     ylabel(['V (cm/s)'])
        %     pause(1)
        %     xlim([-0.5 0.5])
        ylim(ywin)
        datetick('x')
        %xlim([datenum(2016,1,5), max(X)])
        xlim([min(X(isfinite(X))), min(now,max(X(isfinite(X))))])
        %     grid on
        box on
        xlabel('Date')
        ylabel(tit)
    end
end
thscr2png(['Microcat_' stn '_' datestr(now,'yyyy_mm')],'150',savepath)

%% plot Microcats TS
figure(10); clf
clear D
cmap = jet(6);


for n = 1:6
    X = [];
    P = [];
    T = [];
    C = [];
    for i = 1:numel(data)
        hold on
        x = cell2mat(data(i).Microcats(n).Timestamp);
        
        p = data(i).Microcats(n).Pressure;
        t = data(i).Microcats(n).Temperature;
        c = data(i).Microcats(n).Conductivity;
        
        abs(datenum(2016,1,1)-x)<500;
        ii = find(isfinite((x)) & abs(datenum(2016,1,1)-x)<500 & p > 0 & p < 10000);
        X=[X; x(ii)];
        P=[P; p(ii)];
        T=[T; t(ii)];
        C=[C; c(ii)];
        
        if NumberOfDailyMessages==5
            X(X>datenum(2016,1,21))=datenum(2016,1,21);
        end
    end
    
    S = sw_salt(10*C/sw_c3515,T,P);
    TH = sw_ptmp(S,T,P,0);
    
    %plot(X,S,'-','color',cmap(n,:),'linewidth',2)
    %datetick('x')
    %xlim([min(X), max(X)])
    
    % tit = 'Pressure (dBar)';
    
    %     title(['AD ' num2str(n) ', color = time (days)'])
    %     xlabel(['U (cm/s)'])
    %     ylabel(['V (cm/s)'])
    %     pause(1)
    %     xlim([-0.5 0.5])
    %ylim(ywin)
    %xlim([datenum(2016,1,5), max(X)])
    %     grid on
    
    pl(n) = plot(S,TH,'.','color',cmap(n,:));
    sref = 34.4:0.05:34.65;
    D{n} = num2str(round(median(P)));
    plot(sref,fp_t(sref,ones(size(sref))*median(P)),'-','linewidth',2,'color',0.7*[1 1 1])
    plot(sref,thgade(sref,34.49,-2.45,'temp'),'-b')
    plot(sref,thgade(sref,34.51,-2.5,'temp'),'-r')
    box on
    %xlabel('Date')
    %ylabel(tit)
    
    
end

xlim([34.4671, 34.6328]);
ylim([-2.5146,-2.2248]);
denscont(0);
xlabel('Salinity')
ylabel('Potential temperature')
legend(pl,D,'location','northwest')
title(stn)
thscr2png(['TS_' stn '_' datestr(now,'yyyy_mm')],'150',savepath)

%% plot Aquadopps
figure(2); clf
cmap = jet(numel(data));
for n = 1:4
    subplot(2,2,n)
    for i = 1:numel(data)
        hold on
        plot(data(i).Aquadopps(n).U,data(i).Aquadopps(n).V,'.','color',cmap(i,:))
    end
    colorbar
    title(['AD ' num2str(n) ', color = time (days)'])
    xlabel(['U (cm/s)'])
    ylabel(['V (cm/s)'])
    pause(1)
    xlim([-0.5 0.5])
    ylim([-0.5 0.5])
    axis square
    grid on
    box on
end
thscr2png(['Aquadopp_' stn '_' datestr(now,'yyyy_mm')],'150',savepath)

%% plot Aquadopps series
figure(20); clf
cmap = jet(numel(data));
for n = 1:4
    subplot(4,1,n)
    for i = 1:numel(data)
        hold on
        if isempty( data(i).Aquadopps(n).Timestamp{1})
            data(i).Aquadopps(n).Timestamp{1} =nan;
        end
        
        plot(data(i).Aquadopps(n).Timestamp{1}+[0:2:22]/24,sqrt(data(i).Aquadopps(n).U.^2+data(i).Aquadopps(n).V.^2),'.','color',cmap(i,:))
    end
    %     colorbar
    %     title(['AD ' num2str(n) ', color = time (days)'])
    %     xlabel(['U (cm/s)'])
    %     ylabel(['V (cm/s)'])
    %     pause(1)
    %     xlim([-0.5 0.5])
    %     ylim([-0.5 0.5])
    %     axis square
    %     grid on
    %     box on
end
thscr2png(['Aquadopp_ser' stn '_' datestr(now,'yyyy_mm')],'150',savepath)
%%