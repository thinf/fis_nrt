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
        while i <= NumberOfDailyMessages +2
           % disp(i)
           
            % keyboard
            Filenames{i} = strcat(path, sprintf('%06d', (MessageNo + i-1)),filetype);
            FileIds{i} = fopen(Filenames{i}, 'r');
            
            %disp(Filenames{i})
            if FileIds{i} ~= -1
              %  disp('filedids exists')
                TempMsgNo = fread(FileIds{i}, 1, '*uint8', 2);
                TempDate{i} = datestr(floor(double((fread(FileIds{i}, 1, '*uint32'))/24/3600 + datenum(1990,1,1,0,0,0))));
                disp(TempDate{i})
            else
              %  disp('do other break')
                
                i = i + 1;
                
                
                break
            end
            
            %First time round loop stores filename and date in relevent cells in
            %array
            if i == 1;
              %  disp('first one')
                MessageArray{TempMsgNo, 1, DayNo} = TempDate{i};
                MessageArray{TempMsgNo, 2, DayNo} = Filenames{i};
                MessageArray{TempMsgNo, 3, DayNo} = FileIds{i};
                
                %Next times round loop, if date matches, stores filename and date in
                %relevent cells in array
            elseif TempDate{i} == TempDate{i-1}
               % disp('same day')
                MessageArray{TempMsgNo, 1, DayNo} = TempDate{i};
                MessageArray{TempMsgNo, 2, DayNo} = Filenames{i};
                MessageArray{TempMsgNo, 3, DayNo} = FileIds{i};
                
                %As soon as there is a message with new date, we break from this
                %loop.
            else
               % disp('else')
                iMax = i;    % max value of i
                breakflag = true;
                break;
            end
            
            if breakflag == true
               % disp('break')
                break
            end
            
           % disp('shift forward')
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
               % MessageArray{j, 2, DayNo} = strcat(workpath, 'blankmsg', sprintf('%d', j), '.sbd');
               % MessageArray{j, 2, DayNo} = strcat('./', 'blankmsg', sprintf('%d', j), '.sbd');
                MessageArray{j, 2, DayNo} = strcat(blankmsgs, sprintf('%d', j), '.sbd'); % repalce missing message by dummy
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
        
        for m = 1:NumberOfDailyMessages
            fids{m} = MessageArray{m, 3, DayNo};
        end
        