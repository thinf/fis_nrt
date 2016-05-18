function analysis_filchner_nrtfun_plot(path,stn,workpath,savepath)
% First attemtp to read and plot the FISP Iridium data
% still heavily under construction
% thinf, 04.02.2015, tore.hattermann@awi.de
% thinf, 14.04.2016. extended to run as standalone
% thinf, 17.05.2016 add ploting for infoboard, input needs to conatin:

if 0
   path = '/home/csys/thatterm/Dropbox/Osci/FISP/#DATA/SBD/300234061032780_'
   stn = 'fse2'
   workpath = '/home/csys/thatterm/fis_nrt' 
   savepath = '/home/csys/thatterm/fis_nrt/'
end

if 0
path = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061032780_' % path and
% prefix of messages
stn = 'fse2' % name of station, on ly used for saving prefix (I hope)
workpath = 'C:\THINF\#TEMP\fis_nrt_tmp' % Non-dropbox synched working
% directory that also conatins the dummy messages
savepath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\infoboard/'
% path for saving output
end


if 0
path = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061031800_' % path and
% prefix of messages
stn = 'fsw1' % name of station, on ly used for saving prefix (I hope)
workpath = 'C:\THINF\#TEMP\fis_nrt_tmp' % Non-dropbox synched working
% directory that also conatins the dummy messages
savepath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\infoboard/'
% path for saving output
end

%
% Janik:
% hier die beiden Pf�de:
%  
% /hs/datex/ingest/mooring/fispocean/fsw1/data/
% /hs/datex/ingest/mooring/fispocean/fse2/data/
%  
% Es sind nun 2 XML Dateien daf�r zust�ndig. 
% 
% Das Skript liegt aktuell /hs/datex/ingest/mooring/fispocean/ dort. 
% 
%
%
% Change this path to the folder and filename prefix dthat the SBD messages
% are stored in.
% Please note that datalogger at FSW1 SBD messages begin with 300234061031800_
% Please note that datalogger at FSE2 SBD messages begin with 300234061038780_
%
% path = 'E:\SBDs\from Gerd\300234061031800_';
% path = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\800/300234061031800_';
%paths{1} = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061031800_';  stns{1} = 'FSW1'
%paths{2} = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061032780_'; stns{2} = 'FSE2';

%paths{1} = 'data/300234061031800_';  stns{1} = 'FSW1'
%paths{2} = 'data/300234061032780_'; stns{2} = 'FSE2';

% for srcipt only
%for ipi = 1:numel(paths)
%    clearvars('-except','paths','ipi', 'stns')
%    path = paths{ipi}
%    stn = stns{ipi}

% specify folder where temporary files will be created. Also the files
    % "blankmsg1.sbd" to "blankmsg7.sbd" MUST be stored here.
    % workpath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\processed/';
    % workpath = './';
    
    % workpath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\nrt_demo\';
    %  workpath = '/hs/datex/ingest/mooring/fispocean/';
    
    nMC = 6; % number of Microcats
    nAD = 4; % number of Aquadopps
    
    NumberOfDailyMessages = 7;
    % The Filchner dataloggers have 7 SBD messages sent by Iridium every day.
    % Earlier SBDs (ca. before 21.01.2016) from IMEI ...800 had 5 messages only.
    % The program can read these messages as well, but the Aquadopp data from this
    % period will make no sense and the Microcat data will have gaps.
    % To properly read those older files, you can set the value back to 5, but
    % then the Program will fail to read newer files
    
    %######## END OF USER DEFINED SECTION ###############
    %%
    % SBD file extension
    filetype = '.sbd';
    
    % %MaxNumberOfDays - this number should be changed to be equa l to the number
    % %of days data you want to process.  It might be worth changing it to (say)
    % %365 days now, although this will mean the script takes a long time to
    % %run.
    % now determined automatically:
    MaxNumberOfDays = round(now-datenum(2015,12,24));%100;
    files = dir([path '*']);
    
    if NumberOfDailyMessages == 7;
        add_dummybytes = 1; % set to 0 for old AD-format/
    elseif NumberOfDailyMessages == 5;
        add_dummybytes = 0; % set to 0 for old AD-format/
    else
        error('STOP, number of daily messages suggest unknown dataformat')
    end
    
    %set up 'message array' to store date, filename and fileID
    MessageArray = cell(NumberOfDailyMessages, 3, MaxNumberOfDays);
    
    % set up output variable that will contain the actual data
    data = struct([]);
    
    %% loop to find first valid file in directory.  This means that any files
    %that were sent by the modem and then deleted before deployment (eg test
    %messages etc) are ignored.
    for MessageNo=0:(MaxNumberOfDays*NumberOfDailyMessages)
        
        testFilename = strcat(path, sprintf('%06d', MessageNo),filetype);
        
        if fopen(testFilename) ~= -1;
            disp('start with message')
            disp(testFilename)
            break
        end
        
    end
    
    DayNo = 1;
    %keyboard
    %Loop through messages and process data.
    
    while MessageNo < (MaxNumberOfDays * NumberOfDailyMessages) && DayNo <= MaxNumberOfDays;
        
        %initialises arrays each iteration around loop (to ensure they are
        %empty)
        Filenames = cell(1);
        FileIds = cell(1);
        TempMsgNo = cell(1);
        TempDate = cell(1);
        
        %Opens all sequential messages with the same date and puts msg no/date
        %into arrays.
        i = 1;
        breakflag = false;
        while i <= NumberOfDailyMessages
            % keyboard
            Filenames{i} = strcat(path, sprintf('%06d', (MessageNo + i-1)),filetype);
            FileIds{i} = fopen(Filenames{i}, 'r');
            
            if FileIds{i} ~= -1
                TempMsgNo = fread(FileIds{i}, 1, '*uint8', 2);
                TempDate{i} = datestr(floor(double((fread(FileIds{i}, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0))));
                disp(TempDate{i})
            else break
            end
            
            %First time round loop stores filename and date in relevent cells in
            %array
            if i == 1;
                MessageArray{TempMsgNo, 1, DayNo} = TempDate{i};
                MessageArray{TempMsgNo, 2, DayNo} = Filenames{i};
                MessageArray{TempMsgNo, 3, DayNo} = FileIds{i};
                
                %Next times round loop, if date matches, stores filename and date in
                %relevent cells in array
            elseif TempDate{i} == TempDate{i-1}
                MessageArray{TempMsgNo, 1, DayNo} = TempDate{i};
                MessageArray{TempMsgNo, 2, DayNo} = Filenames{i};
                MessageArray{TempMsgNo, 3, DayNo} = FileIds{i};
                
                %As soon as there is a message with new date, we break from this
                %loop.
            else
                iMax = i;    % max value of i
                breakflag = true;
                break;
            end
            
            if breakflag == true
                break
            end
            
            iMax = i;
            i = i + 1;
            
        end
        
        %records MessageNo for next iteration round loop.
        MessageNo = MessageNo+i-1;
        %MessageNo = MessageNo+i;
        
        %resets each file in turn.  If any one is missing, replaces with a blank
        %file.  Missing messages can be caused by the Iridium Constellation not
        %receiving messages from the datalogger correctly.
        for j = 1:1:NumberOfDailyMessages;
            
            if isempty(MessageArray{j, 3, DayNo})
                %MessageArray{j, 2, DayNo} = strcat('C:\Users\THOSTR\Documents\documents from field season 2015-16\inductive system\FSW1\Iridium Analysis\blankmsg', sprintf('%d', j), '.sbd');
                MessageArray{j, 2, DayNo} = strcat(workpath, '/blankmsg', sprintf('%d', j), '.sbd');
              %  MessageArray{j, 2, DayNo} = strcat('./', 'blankmsg', sprintf('%d', j), '.sbd');
                MessageArray{j, 3, DayNo} = fopen(MessageArray{j, 2, DayNo}, 'r');
                
            else
                
                %check file size - if file is less than 337 bytes, append file
                %with extra zeros.  This ensures that the daily_sbd function
                %always manages data with the correct number of bytes of data,
                %even if some of the data is 'zero' fields.
                fseek(MessageArray{j, 3, DayNo}, 0 , 'eof');
                position = ftell(MessageArray{j, 3, DayNo});
                
                if position ~= 337;
                    
                    fclose(MessageArray{j, 3, DayNo});
                    MessageArray{j, 3, DayNo} = fopen(MessageArray{j, 2, DayNo}, 'a+');
                    
                    noBytes = 337 - position;
                    
                    for k = 1:1:noBytes;
                        fwrite(MessageArray{j, 3, DayNo}, zeros(1), 'uint8');
                    end
                    
                end
                
                fseek(MessageArray{j, 3, DayNo}, 0 , 'eof');
                position = ftell(MessageArray{j, 3, DayNo});
                
                fclose(MessageArray{j, 3, DayNo});
                MessageArray{j, 3, DayNo} = fopen(MessageArray{j, 2, DayNo}, 'r');
                
            end
            
        end
        
        
        %Strip data from files and put in 'data' structure.
        %     [data(DayNo).Housekeeping, data(DayNo).Microcat1, data(DayNo).Microcat2,...
        %         data(DayNo).Microcat3, data(DayNo).Microcat4, data(DayNo).Microcat5,...
        %         data(DayNo).Microcat6, data(DayNo).Aquadopp1, data(DayNo).Aquadopp2,...
        %         data(DayNo).Aquadopp3, data(DayNo).Aquadopp4]...
        %         = dailySBD_filchner_FSW1_th(MessageArray{1, 3, DayNo},...
        %         MessageArray{2, 3, DayNo}, MessageArray{3, 3, DayNo},...
        %         MessageArray{4, 3, DayNo}, MessageArray{5, 3, DayNo}, workpath);
        %
        
        for m = 1:NumberOfDailyMessages
            fids{m} = MessageArray{m, 3, DayNo};
        end
        % [data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
        %         = dailySBD_filchner_FSW1_th(nMC, nAD, MessageArray{1, 3, DayNo},...
        %         MessageArray{2, 3, DayNo}, MessageArray{3, 3, DayNo},...
        %         MessageArray{4, 3, DayNo}, MessageArray{5, 3, DayNo}, workpath);
        
        [data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
            = dailySBD_filchner(nMC, nAD, fids, workpath, add_dummybytes);
        
        %% write out data of the day
        if 0
        if isfinite(data(DayNo).Housekeeping.Date)
           % dfile = [paths{ipi} stns{ipi} '_' datestr(datenum(data(DayNo).Housekeeping.Date),'yyyy-mm-dd')];
            dfile = [savepath stn '_' datestr(datenum(data(DayNo).Housekeeping.Date),'yyyy-mm-dd')];
            %%
            if ~exist([dfile '.hdr'],'file')
                diary([dfile '.hdr'])
                disp(data(DayNo).Housekeeping)
                diary off
                
                
                %%  Write Microcat data
                fields = {'Timestamp','Pressure','Temperature','Conductivity'};
                units = {'JulianDays','dBar','degC','mS/cm'};
                suffix = '_Microcats';
                file = [dfile suffix '.txt'];
                fid = fopen(file,'w+');
                dmat = [];
                for mcs = 1:nMC;
                    % write header
                    
                    for iif = 1:numel(fields)
                        fprintf(fid,'%s',['MC' num2str(mcs) ':' fields{iif} ':' units{iif} ';']);
                    end
                    dmat = [dmat, [data(DayNo).Microcats(mcs).Timestamp{:}]',...
                        data(DayNo).Microcats(mcs).Pressure, ...
                        data(DayNo).Microcats(mcs).Temperature, ...
                        data(DayNo).Microcats(mcs).Conductivity*10,...
                        ];
                    
                end
                
                fprintf(fid,'\n');
                fclose(fid);
                
                dlmwrite(file,dmat,'delimiter',';','precision',9,'-append')
                
                %%  Write Aquadopp data
                
                
                fields = {'Timestamp','U','V','W','Head','Pitch','Roll','Pressure','Temperature'};
                units = {'JulianDays','m/s','m/s','m/s','deg','deg','deg','dBar','degC'};
                suffix = '_Aquadopps';
                file = [dfile suffix '.txt'];
                fid = fopen(file,'w+');
                dmat = [];
                
                for mcs = 1:nAD;
                    % write header
                    
                    for iif = 1:numel(fields)
                        fprintf(fid,'%s',['AD' num2str(mcs) ':' fields{iif} ':' units{iif} ';']);
                    end
                    
                    adtime = [data(DayNo).Aquadopps(mcs).Timestamp{:}]';
                    
                    %if isempty( data(DayNo).Aquadopps(mcs).Timestamp{1})
                    %   data(DayNo).Aquadopps(mcs).Timestamp =-9999;
                    %end
                    %num_ = data(i).Aquadopps(n).Timestamp{1}+[0:2:22]'/24; %cell2mat(data(i).Aquadopps(n).Timestamp);
                    adtime(isnan(adtime))=-9999;
                    
                    dmat = [dmat, adtime,...
                        data(DayNo).Aquadopps(mcs).U, ...
                        data(DayNo).Aquadopps(mcs).V, ...
                        data(DayNo).Aquadopps(mcs).W, ...
                        data(DayNo).Aquadopps(mcs).Head, ...
                        data(DayNo).Aquadopps(mcs).Pitch, ...
                        data(DayNo).Aquadopps(mcs).Roll, ...
                        data(DayNo).Aquadopps(mcs).P, ...
                        data(DayNo).Aquadopps(mcs).T, ...
                        ];
                end
                
                fprintf(fid,'\n');
                fclose(fid);
                dlmwrite(file,dmat,'delimiter',';','precision',9,'-append')
            end   
        end
        end % if 0 writ out data of the day
        %%
        DayNo = DayNo + 1;
    end
% end % run over paths for script only

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
            ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<500 & p_ > 0 & p_ < 10000);
            
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
            ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<500 & abs(cv_)<1);
            
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
numu = [datenum(2016,01,06):1:round(now)-1];
clear tu pu

if strcmp(stn(1:4),'fsw1')
    ni = 1:numel(mc.num)-1;
else
    ni = 1:numel(mc.num);
end

for i = ni
    ii = i + numel(mc.num)-numel(ni);
    tu(:,i) = interp1(mc.num{ii},vfilt(mc.t{ii},12),numu);
    pu(i) = median(mc.p{ii});
end
sp1 = subplot(4,1,1:2);
contourf(numu,pu,tu',24,'color','none'); shading flat;
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



sp2=subplot(4,1,3:4);
sp1pos =  get(sp1,'position');
sp2pos =  get(sp2,'position');
sp2pos(3) = sp1pos(3)*0.925;

%
 adp =[6 5 2 1];
 clear adpu
 for i = 1:numel(ad.num)
    adpu{i} = [num2str(round(median(mc.p{adp(i)}))) ' dBar'];
end

    p1 = plot(ad.num{1},abs(ad.cv{1})*100,'.',...
        ad.num{2},abs(ad.cv{2})*100,'.',...
        ad.num{3},abs(ad.cv{3})*100,'.',...
        ad.num{4},abs(ad.cv{4})*100,'.');
    
    hold on
    fl = 60;
    p2 = plot(ad.num{1},vfilt(abs(ad.cv{1})*100,fl),...
        ad.num{2},vfilt(abs(ad.cv{2})*100,fl),...
        ad.num{3},vfilt(abs(ad.cv{3})*100,fl),...
        ad.num{4},vfilt(abs(ad.cv{4})*100,fl),'linewidth',3);
    set(gca,'ylim',[0 50])
    ylabel('Speed [cm/s]','fontsize',16)
    set(sp2,'xlim',get(sp1,'xlim'),'fontsize',16)
    grid on
    datetick('x','keeplimits')
    lgd = legend(p2,adpu,'location','eastoutside');
set(sp2,'position',sp2pos)

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
exit
end
