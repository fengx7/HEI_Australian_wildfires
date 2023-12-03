%% Read data from the original files
for iyr  = 1:14
    year_num = 200800 + (iyr-1)*100;
for imon = 1:12
    date_tmp = num2str(year_num+imon);
    filename = ['E:\Work_for_wildfires\MODIS_C6\MCD14ML\MCD14ML.',date_tmp,'.006.03.txt'];
    delimiter = ' ';
    startRow = 2;

    % Format:
    formatSpec = '%f%f%s%f%f%f%f%f%f%f%f%s%[^\n\r]';

    % Open file
    fileID = fopen(filename,'r');

    % Read data
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    % Close file
    fclose(fileID);

    % Save variables
    YYYYMMDD = dataArray{:, 1};
    HHMM = dataArray{:, 2};
    sat = dataArray{:, 3};
    lat = dataArray{:, 4};
    lon = dataArray{:, 5};
    T21 = dataArray{:, 6};
    T31 = dataArray{:, 7};
    sample = dataArray{:, 8};
    FRP = dataArray{:, 9};
    conf = dataArray{:, 10};
    type1 = dataArray{:, 11};
    dn = dataArray{:, 12};
    save_filename = [filename(1:end-3),'mat'];
    disp(save_filename);
    save(save_filename,'YYYYMMDD','HHMM','sat','lat','lon','T21','T31','sample','FRP','conf','type1','dn');
    
    % clear variables
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
end
end



