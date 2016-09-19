 % %MaxNumberOfDays - this number should be changed to be equal to the number
    % %of days data you want to process.  It might be worth changing it to (say)
    % %365 days now, although this will mean the script takes a long time to
    % %run.
    % now determined automatically:
    MaxNumberOfDays = round(now-datenum(2015,12,24));%100;
    files = dir([path '*']);
    
    %set up 'message array' to store date, filename and fileID
    MessageArray = cell(NumberOfDailyMessages, 3, MaxNumberOfDays);
    
    % initialize output variable that will contain the actual data
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