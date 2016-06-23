function [Housekeeping, Microcats, Aquadopps]...
    = dailySBD_filchner(nMC, nAD, fids, workpath, msg_type)

% % add reference date for checking
% if exist([workpath 'prevday'],'file')
%     load([workpath 'prevday'], 'LoggerTimeNum')
%     oldLoggerTimeNum = LoggerTimeNum;
%     dtmax = 5;
% else
%     oldLoggerTimeNum = 0.5*datenum(2015,12,24)+now;
%     dtmax = 0;
% end
dtmax = abs(datenum(2015,12,24)-now); % time stamp threshold

%generate data structures for storing results
Housekeeping = struct('Date', [], 'MsgNo', [], 'MsgTimestamps', [],...
    'BattV',[],'Temp',[], 'Latitude', [], 'Longitude',[],'VM2Time', [], 'GPSTime',...
    [],'TotalSamples', [], 'SDCount', [], 'ResetSource', [], 'AquadoppsBattV',[]);

for i = 1:nAD
    Aquadopps(i) = struct('Timestamp', [], 'U', [], 'V', [], 'W', [], 'LoggerTime',[]);
end

for i = 1:nMC
    Microcats(i) = struct('Timestamp', [], 'Conductivity', [], 'Temperature', [], 'Pressure', [], 'LoggerTime',[]);
end
%Save Measurment Interval in seconds here.  The current deployment has
%the instruments running every 2 hours (7200 seconds).
MeasurementInt = 7200;

%discard first byte - we have finished with this.  This byte contains
%the message number which we used in 'analysis_petermann1.m' to
%identify the messages.
%  keyboard
for m = 1:numel(fids)
    % fids{m}
    fread(fids{m}, 1, '*uint8');
    Housekeeping.MsgNo(m) = fread(fids{m}, 1, '*uint16');
    MsgTimestamps{m} = datestr(double(fread(fids{m}, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0));
end
Housekeeping.Date =  datestr(floor(datenum(datevec(MsgTimestamps{1}))));
% fread(fid2, 1, '*uint8');
% fread(fid3, 1, '*uint8');
% fread(fid4, 1, '*uint8');
% fread(fid5, 1, '*uint8');

%read message numbers

% Housekeeping.MsgNo2 = fread(fid2, 1, '*uint16');
% Housekeeping.MsgNo3 = fread(fid3, 1, '*uint16');
% Housekeeping.MsgNo4 = fread(fid4, 1, '*uint16');
% Housekeeping.MsgNo5 = fread(fid5, 1, '*uint16');

%read timestamps
% MsgTimestamp1 = double(fread(fid1, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0);
% MsgTimestamp2 = double(fread(fid2, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0);
% MsgTimestamp3 = double(fread(fid3, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0);
% MsgTimestamp4 = double(fread(fid4, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0);
% MsgTimestamp5 = double(fread(fid5, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0);

%convert timestamps to datestrings.
% Housekeeping.Date = datestr(floor(MsgTimestamp1));
% Housekeeping.MsgTimestamp1 = datestr(MsgTimestamp1);
% Housekeeping.MsgTimestamp2 = datestr(MsgTimestamp2);
% Housekeeping.MsgTimestamp3 = datestr(MsgTimestamp3);
% Housekeeping.MsgTimestamp4 = datestr(MsgTimestamp4);
% Housekeeping.MsgTimestamp5 = datestr(MsgTimestamp5);

%Recreate original binary file from remainder of messages (headers now
%removed). This original binary file is the daily data buffer which is
%stored in the memory of the datalogger before being sent by Iridium in
%337 byte segments.
BinaryArray = [];
for m = 1:numel(fids)
    if m==1
        len = 500;
    else
        len = 340;
    end
    BinaryArray = cat(1, BinaryArray, fread(fids{m}, len, '*uint8'));
end
%BinaryArray = cat(1, BinaryArray, fread(fids{m}, 500, '*uint8'), fread(fid2, 340, '*uint8'),...
%    fread(fid3, 340, '*uint8'), fread(fid4, 340, '*uint8'), fread(fid5, 340, '*uint8'));

% TempFileID = fopen('C:\Users\THOSTR\Documents\MATLAB\tempbinarray.bin','w');
TempFileID = fopen([workpath '\tempbinarray.bin'],'w');
fwrite(TempFileID, BinaryArray);
%TempFileID = fopen('C:\Users\THOSTR\Documents\MATLAB\tempbinarray.bin', 'r');
TempFileID = fopen([workpath '\tempbinarray.bin'], 'r');

%Get housekeeping, location and time data
Housekeeping.BattV = double(fread(TempFileID, 1, '*int16'))/100;

%note that this temp field is now 32 bit - used to be 16 bit on
%previous (Ronne) systems.
Housekeeping.Temp = double(fread(TempFileID, 1, '*int32'))/100;

Housekeeping.Latitude = fread(TempFileID, 1, '*int32');
tmp = dec2hex(Housekeeping.Latitude);
if length(tmp) > 5
    Housekeeping.Latitude = double(hex2dec(tmp(1:5)))/10000;
else
    Housekeeping.Latitude = double(hex2dec(tmp))/10000;
end

Housekeeping.Longitude = fread(TempFileID, 1, '*int32');
tmp = dec2hex(Housekeeping.Longitude);
if length(tmp) > 5
    Housekeeping.Longitude = double(hex2dec(tmp(1:5)))/10000;
else
    Housekeeping.Longitude = double(hex2dec(tmp))/10000;
end
%keyboard
%Housekeeping.Longitude = double(fread(TempFileID, 1, '*int32'))/10000;

Housekeeping.GPSTime = double(fread(TempFileID, 1, '*int32'));
Housekeeping.VM2Time = double(fread(TempFileID, 1, '*int32'));

Housekeeping.TotalSamples = swapbytes(fread(TempFileID, 1, '*int16'));
Housekeeping.SDCount = fread(TempFileID, 1, '*int8');
Housekeeping.ResetSource = fread(TempFileID, 1, '*int8');

for i = 1:nAD
    Housekeeping.AquadoppsBattV(i) = double(fread(TempFileID, 1, '*int16'))/10;
end

%replace Zeros with NaNs (DateNo 726834 refers to 01-Jan-1990 epoch).
if datenum(Housekeeping.Date) == 726834, Housekeeping.Date = NaN; end
for m = 1:numel(fids)
    if Housekeeping.MsgNo(m) == 0, Housekeeping.MsgNo(m) = NaN; end
    %  if datenum(Housekeeping.MsgTimestamps(m,:)) == 726834, Housekeeping.MsgTimestamps(m,:) = NaN; end
end
% if Housekeeping.MsgNo2 == 0, Housekeeping.MsgNo2 = NaN; end
% if datenum(Housekeeping.MsgTimestamp2) == 726834, Housekeeping.MsgTimestamp2 = NaN; end
% if Housekeeping.MsgNo3 == 0, Housekeeping.MsgNo3 = NaN; end
% if datenum(Housekeeping.MsgTimestamp3) == 726834, Housekeeping.MsgTimestamp3 = NaN; end
% if Housekeeping.MsgNo4 == 0, Housekeeping.MsgNo4 = NaN; end
% if datenum(Housekeeping.MsgTimestamp4) == 726834, Housekeeping.MsgTimestamp4 = NaN; end
% if Housekeeping.MsgNo5 == 0, Housekeeping.MsgNo4 = NaN; end
% if datenum(Housekeeping.MsgTimestamp5) == 726834, Housekeeping.MsgTimestamp4 = NaN; end
if Housekeeping.BattV == 0, Housekeeping.BattV = NaN; end
if Housekeeping.Temp == 0, Housekeeping.Temp = NaN; end
if Housekeeping.Latitude == 0, Housekeeping.Latitude = NaN; end
if Housekeeping.Longitude == 0, Housekeeping.Longitude = NaN; end
if Housekeeping.GPSTime == 0, Housekeeping.GPSTime = NaN; end
if Housekeeping.VM2Time == 0, Housekeeping.VM2Time = NaN; end
if Housekeeping.TotalSamples == 0, Housekeeping.TotalSamples = NaN; end
if Housekeeping.SDCount == 0, Housekeeping.SDCount = NaN; end
if Housekeeping.ResetSource == 0, Housekeeping.ResetSource = NaN; end
for i = 1:nAD
    if Housekeeping.AquadoppsBattV(i) == 0, Housekeeping.AquadoppsBattV(i) = NaN; end
end

%Cycle through remianing data to get data from the instruments.
Repeat = (24*3600)/MeasurementInt;
LoggerTimeNum = -9999;
for j= 1:1:Repeat
    
    %**TO DO!  Chekc if we are in first 2hours of day and runsome
    %checks on the final timestamp.  We could be using the
    %dtaloggertime to get the date from one day, with the sample time
    %in the previous day.
    
    % currently: fix day from logger timestamp and add sample time from
    % insturments. Also, originally this is done using the datetime
    % function, which has been introduced after the date of Tore's
    % matlab license, thus he had to modify it to use the old datenum
    % format.
    
    Time = double(fread(TempFileID, 1, '*int32'));
    
    % thinf: use datenum/ datevec instead
    % LoggerTimeStamp = datetime(datestr((Time)/24/3600 + datenum(1990,1,1,0,0,0)));
    % Year = year(LoggerTimeStamp);
    % Month = month(LoggerTimeStamp);
    % Day = day(LoggerTimeStamp);
    
    oldLoggerTimeNum = LoggerTimeNum;
    LoggerTimeNum = (Time)/24/3600 + datenum(1990,1,1,0,0,0);
    
    if oldLoggerTimeNum~=-9999 && abs(LoggerTimeNum-oldLoggerTimeNum) > 3
        LoggerTimeNum = -9999; %oldLoggerTimeNum;
    elseif abs(LoggerTimeNum-now)> dtmax
        LoggerTimeNum = -9999;
    end
    
    %LoggerTimeVec = datevec((Time)/24/3600 + datenum(1990,1,1,0,0,0));
    LoggerTimeVec = datevec(LoggerTimeNum);
    Year = LoggerTimeVec(1);
    Month = LoggerTimeVec(2);
    Day = LoggerTimeVec(3);
    
    
    for i = 1:nMC
        Hour = double(fread(TempFileID, 1, '*int8'));
        Min = double(fread(TempFileID, 1, '*int8'));
        Sec = double(fread(TempFileID, 1, '*int8'));
        
        mctimenum = datenum(Year, Month, Day, Hour, Min, Sec);
         if ~isfinite(mctimenum) || abs(datenum(2016,1, 1)-mctimenum)> dtmax
            mctimenum = -9999;
         end
        
        Microcats(i).Timestamp{j,1} = mctimenum;
        % assume that we are reading data from previous day if MC timestamp
        % is younger than logger timestamp (not necessairily true for large
        % clock drift)
        
        if Microcats(i).Timestamp{j,1} > LoggerTimeNum
            Microcats(i).Timestamp{j,1} = Microcats(i).Timestamp{j,1}-1;
        end
        Microcats(i).Conductivity(j, 1) = double(swapbytes(fread(TempFileID, 1, '*int32')))/100000;
        Microcats(i).Temperature(j, 1) = -1 * double(swapbytes(fread(TempFileID, 1, '*int16')))/10000;
        Microcats(i).Pressure(j, 1) = double(swapbytes(fread(TempFileID, 1, '*int32')))/1000;
        Microcats(i).LoggerTime{j, 1} = LoggerTimeNum;
    end
    
    for i = 1:nAD
        %        if j==10 && i ==2; keyboard; end
        % bdc2deg: defined below, workaround is str2num(char(dec2hex(bdc})))
        % Hour = double(bcd2dec(fread(TempFileID, 1, '*int8')));
        % Min = double(bcd2dec(fread(TempFileID, 1, '*int8')));
        % Sec = double(bcd2dec(fread(TempFileID, 1, '*int8')));
        
        Hour =nan;
        Min =nan;
        Sec = nan;
 %       keyboard
        try
            Hour = str2num(char(dec2hex(fread(TempFileID, 1, '*int8'))));
            Min = str2num(char(dec2hex(fread(TempFileID, 1, '*int8'))));
            Sec = str2num(char(dec2hex(fread(TempFileID, 1, '*int8'))));
        
        if Hour>24 || Hour < 0
            Hour=[];
            Min =[];
            Sec=[];
        end
        end
        
        clear adtimedum
        Min(isempty(Min))=0;
        Sec(isempty(Sec))=0;
        Hour(isempty(Hour))=j*2/22;
        adtimedum = datenum(Year, Month, Day, Hour, Min, Sec);
        adtimedum(isempty(adtimedum)) = -9999;
        
        if ~isfinite(adtimedum) || abs(datenum(2016,1, 1)-adtimedum)> dtmax
            adtimedum = -9999;
        end
        
        if adtimedum > LoggerTimeNum;
            adtimedum = adtimedum-1;
        end
        
        Aquadopps(i).Timestamp{j,1} = adtimedum;
        Aquadopps(i).U(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/1000; % eastward, m/s
        Aquadopps(i).V(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/1000;
        Aquadopps(i).W(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/1000;
        
        if msg_type==1
        Aquadopps(i).P(j, 1) = double(int16(65536*fread(TempFileID, 1, '*int8')) + (fread(TempFileID, 1, '*int16', 0, 'l')))/1000; % dBar
        Aquadopps(i).T(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/100; % degC
        Aquadopps(i).Head(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/10; %deg
        Aquadopps(i).Pitch(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/10; %deg
        Aquadopps(i).Roll(j, 1) = double(fread(TempFileID, 1, '*int16', 0, 'l'))/10; %deg
        Aquadopps(i).Amp1(j, 1) = (fread(TempFileID, 1, '*int8')); %counts
        Aquadopps(i).Amp2(j, 1) = (fread(TempFileID, 1, '*int8'));
        Aquadopps(i).Amp3(j, 1) = (fread(TempFileID, 1, '*int8'));
        end
        Aquadopps(i).LoggerTime{j,1} = LoggerTimeNum;
        
%         if add_dummybytes==1
%             % disp('long AD format, add dummybytes')
%             % read up to 20 bytes AD data
%             dum1 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%             dum2 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%             dum3 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%             dum4 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%             dum5 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%             dum6 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%             dum7 = double(fread(TempFileID, 1, '*int16', 0, 'l'));
%         end
        
    end
    
    %         %Hour = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         %Min = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         %Sec = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Hour = fread(TempFileID, 1, '*int8');
    %         Min = fread(TempFileID, 1, '*int8');
    %         Sec = fread(TempFileID, 1, '*int8');
    %         %Aquadopp2.Timestamp{j,1} = datetime(Year, Month, Day, Hour, Min, Sec);
    %         Aquadopp2.Timestamp{j,1} = Hour;
    %         %if Aquadopp2.Timestamp{j,1} > LoggerTimeStamp, Aquadopp2.Timestamp{j,1} = Aquadopp2.Timestamp{j,1}-datenum(0000,0,1,0,0,0); end
    %         Aquadopp2.XCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp2.YCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp2.ZCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp2.LoggerTime{j,1} = LoggerTimeStamp;
    %
    %         Hour = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Min = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Sec = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Aquadopp3.Timestamp{j,1} = datetime(Year, Month, Day, Hour, Min, Sec);
    %          if Aquadopp3.Timestamp{j,1} > LoggerTimeStamp, Aquadopp3.Timestamp{j,1} = Aquadopp3.Timestamp{j,1}-datenum(0000,0,1,0,0,0); end
    %         Aquadopp3.XCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp3.YCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp3.ZCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp3.LoggerTime{j,1} = LoggerTimeStamp;
    %
    %         Hour = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Min = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Sec = bcd2dec(fread(TempFileID, 1, '*int8'));
    %         Aquadopp4.Timestamp{j,1} = datetime(Year, Month, Day, Hour, Min, Sec);
    %          if Aquadopp4.Timestamp{j,1} > LoggerTimeStamp, Aquadopp4.Timestamp{j,1} = Aquadopp4.Timestamp{j,1}-datenum(0000,0,1,0,0,0); end
    %         Aquadopp4.XCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp4.YCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp4.ZCurrent(j, 1) = (fread(TempFileID, 1, '*int16', 0, 'l'));
    %         Aquadopp4.LoggerTime{j,1} = LoggerTimeStamp;
    
    %For any zero fields, replace with NaNs.  (DateNo 726834 refers to 01-Jan-1990 epoch).
    for i=1:nMC
        if datenum(Microcats(i).Timestamp{j, 1}) < 726836, Microcats(i).Timestamp{j, 1} = NaN; end
        if Microcats(i).Conductivity(j, 1) == 0, Microcats(i).Conductivity(j, 1) = NaN; end
        if Microcats(i).Temperature(j, 1) == 0, Microcats(i).Temperature(j, 1) = NaN; end
        if Microcats(i).Pressure(j, 1) == 0 || Microcats(i).Pressure(j, 1) == 500, Microcats(i).Pressure(j, 1) = NaN; end
        if datenum(Microcats(i).LoggerTime{j, 1}) < 726836, Microcats(i).LoggerTime{j, 1} = NaN; end
    end
    %         if datenum(Microcat2.Timestamp{j, 1}) < 726836, Microcat2.Timestamp{j, 1} = NaN; end
    %         if Microcat2.Conductivity(j, 1) == 0, Microcat2.Conductivity(j, 1) = NaN; end
    %         if Microcat2.Temperature(j, 1) == 0, Microcat2.Temperature(j, 1) = NaN; end
    %         if Microcat2.Pressure(j, 1) == 0, Microcat2.Pressure(j, 1) = NaN; end
    %         if datenum(Microcat2.LoggerTime{j, 1}) < 726836, Microcat2.LoggerTime{j, 1} = NaN; end
    %         if datenum(Microcat3.Timestamp{j, 1}) < 726836, Microcat3.Timestamp{j, 1} = NaN; end
    %         if Microcat3.Conductivity(j, 1) == 0, Microcat3.Conductivity(j, 1) = NaN; end
    %         if Microcat3.Temperature(j, 1) == 0, Microcat3.Temperature(j, 1) = NaN; end
    %         if Microcat3.Pressure(j, 1) == 0, Microcat3.Pressure(j, 1) = NaN; end
    %         if datenum(Microcat3.LoggerTime{j, 1}) < 726836, Microcat3.LoggerTime{j, 1} = NaN; end
    %         if datenum(Microcat4.Timestamp{j, 1}) < 726836, Microcat4.Timestamp{j, 1} = NaN; end
    %         if Microcat4.Conductivity(j, 1) == 0, Microcat4.Conductivity(j, 1) = NaN; end
    %         if Microcat4.Temperature(j, 1) == 0, Microcat4.Temperature(j, 1) = NaN; end
    %         if Microcat4.Pressure(j, 1) == 0, Microcat4.Pressure(j, 1) = NaN; end
    %         if datenum(Microcat4.LoggerTime{j, 1}) < 726836, Microcat4.LoggerTime{j, 1} = NaN; end
    %         if datenum(Microcat5.Timestamp{j, 1}) < 726836, Microcat5.Timestamp{j, 1} = NaN; end
    %         if Microcat5.Conductivity(j, 1) == 0, Microcat5.Conductivity(j, 1) = NaN; end
    %         if Microcat5.Temperature(j, 1) == 0, Microcat5.Temperature(j, 1) = NaN; end
    %         if Microcat5.Pressure(j, 1) == 0, Microcat5.Pressure(j, 1) = NaN; end
    %         if datenum(Microcat5.LoggerTime{j, 1}) < 726836, Microcat5.LoggerTime{j, 1} = NaN; end
    %         if datenum(Microcat6.Timestamp{j, 1}) < 726836, Microcat6.Timestamp{j, 1} = NaN; end
    %         if Microcat6.Conductivity(j, 1) == 0, Microcat6.Conductivity(j, 1) = NaN; end
    %         if Microcat6.Temperature(j, 1) == 0, Microcat6.Temperature(j, 1) = NaN; end
    %         if Microcat6.Pressure(j, 1) == 0, Microcat6.Pressure(j, 1) = NaN; end
    %         if datenum(Microcat6.LoggerTime{j, 1}) < 726836, Microcat6.LoggerTime{j, 1} = NaN; end
    
    for i=1:nAD
        if datenum(Aquadopps(i).Timestamp{j, 1}) < 726836, Aquadopps(i).Timestamp{j, 1} = NaN; end
        if Aquadopps(i).U(j, 1) == 0, Aquadopps(i).U(j, 1) = NaN; end
        if Aquadopps(i).V(j, 1) == 0, Aquadopps(i).V(j, 1) = NaN; end
        if Aquadopps(i).W(j, 1) == 0, Aquadopps(i).W(j, 1) = NaN; end
        if datenum(Aquadopps(i).LoggerTime{j, 1}) < 726836, Aquadopps(i).LoggerTime{j, 1} = NaN; end
    end
    %         %if datenum(Aquadopp2.Timestamp{j, 1}) < 726836, Aquadopp2.Timestamp{j, 1} = NaN; end
    %         if Aquadopp2.XCurrent(j, 1) == 0, Aquadopp2.XCurrent(j, 1) = NaN; end
    %         if Aquadopp2.YCurrent(j, 1) == 0, Aquadopp2.YCurrent(j, 1) = NaN; end
    %         if Aquadopp2.ZCurrent(j, 1) == 0, Aquadopp2.ZCurrent(j, 1) = NaN; end
    %         if datenum(Aquadopp2.LoggerTime{j, 1}) < 726836, Aquadopp2.LoggerTime{j, 1} = NaN; end
    %         if datenum(Aquadopp3.Timestamp{j, 1}) < 726836, Aquadopp3.Timestamp{j, 1} = NaN; end
    %         if Aquadopp3.XCurrent(j, 1) == 0, Aquadopp3.XCurrent(j, 1) = NaN; end
    %         if Aquadopp3.YCurrent(j, 1) == 0, Aquadopp3.YCurrent(j, 1) = NaN; end
    %         if Aquadopp3.ZCurrent(j, 1) == 0, Aquadopp3.ZCurrent(j, 1) = NaN; end
    %         if datenum(Aquadopp3.LoggerTime{j, 1}) < 726836, Aquadopp3.LoggerTime{j, 1} = NaN; end
    %         if datenum(Aquadopp4.Timestamp{j, 1}) < 726836, Aquadopp4.Timestamp{j, 1} = NaN; end
    %         if Aquadopp4.XCurrent(j, 1) == 0, Aquadopp4.XCurrent(j, 1) = NaN; end
    %         if Aquadopp4.YCurrent(j, 1) == 0, Aquadopp4.YCurrent(j, 1) = NaN; end
    %         if Aquadopp4.ZCurrent(j, 1) == 0, Aquadopp4.ZCurrent(j, 1) = NaN; end
    %         if datenum(Aquadopp4.LoggerTime{j, 1}) < 726836, Aquadopp4.LoggerTime{j, 1} = NaN; end
    %
end

%Close files and delete temp file.
fclose('all');
delete([workpath '\tempbinarray.bin']);

% save([workpath 'prevday'], 'LoggerTimeNum')

end


function [dec] = bcd2dec(bcd)

%hex = dec2hex(bcd);
dec = 10* bitshift(bcd, -4) + bitand(bcd, 16);

end

