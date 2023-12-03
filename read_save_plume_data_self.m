
filepath = 'E:\work_for_wildfires\injecton_height\MISR\';
diroutput = dir(fullfile(filepath,'Plumes*'));
filename = {diroutput.name};

for ii = 1:length(filename)
    filename_tmp = char(filename(ii));
    disp(filename_tmp);
    delimiter = ' ';

    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

    % Open file
    fileID = fopen([filepath,filename_tmp],'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);

    % Close file
    fclose(fileID);

    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
        rawData = dataArray{col};
        for row=1:size(rawData, 1)
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData(row), regexstr, 'names');
                numbers = result.numbers;

                invalidThousandsSeparator = false;
                if numbers.contains(',')
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'))
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end。
                if ~invalidThousandsSeparator
                    numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch
                raw{row, col} = rawData{row};
            end
        end
    end


    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);
    raw(R) = {NaN};

    % Create output variables
    Plumes = cell2mat(raw);
    % Clear variables
    clearvars delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
    
    record_num = Plumes(end,1);
    if ~isnan(record_num)
        Plumes_data = Plumes(end - record_num + 1:end,1:end-2);
        delimiter = ' ';
        formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
        fileID = fopen([filepath,filename_tmp],'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
        fclose(fileID);
        Plumes_s = [dataArray{1:end-1}];

        retrieval_quality = Plumes_s(39,5);
        retrieval_date    = Plumes_s(4,4);
        retrieval_time    = Plumes_s(5,4);
        retrieval_biome   = Plumes_s(25,6:end);
        retrieval_cover_perc = Plumes_s(31,5);
        fp_lon            = Plumes_s(22,5);
        fp_lat            = Plumes_s(23,5);
        fire_elev         = Plumes_s(32,7);
        median_ht         = Plumes_s(33,7);
        max_ht            = Plumes_s(34,7);
        total_frp         = Plumes_s(38,6);
        clearvars delimiter formatSpec fileID dataArray ans;
        variable_name = [filepath,'\MISR_mat\',filename_tmp(1:end-4),'.mat'];
        save(variable_name,'Plumes_data','retrieval_quality','retrieval_biome','retrieval_cover_perc','retrieval_time','retrieval_date','fp_lon','fp_lat','fire_elev','median_ht','max_ht','total_frp');
    end
end
