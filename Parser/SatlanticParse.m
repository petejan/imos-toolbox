function sample_data = SatlanticParse( fn, mode )
%IsusSunaParse Parse a zip file with satlantic ISUS/SUNA data and
%calibration files
%
% Outputs:
%   sample_data - contains a time vector (in matlab numeric format)
%
% Author: 		Peter Jansen
%

%
% Copyright (c) 2017, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS).
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%     * Redistributions of source code must retain the above copyright notice, 
%       this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in the 
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the AODN/IMOS nor the names of its contributors 
%       may be used to endorse or promote products derived from this software 
%       without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%

% Check input, set up data structures
narginchk(1, 2);

if ~iscellstr(fn)
    filenames{1} = fn;
else
    filenames = fn;
end

if ~exist(filenames{1}, 'file')
    e = sprintf('File not found : %s\n', filenames{1});
    error(e);
end

% struct to map satlantic variable names to IMOS variable names

satlanticVar = {'NITRATE','NITRATE_UM', 'AUX1','AUX2','AUX3','RMSe','T_INT','T_SPEC',...
                'T_LAMP','LAMP_TIME','HUMIDITY','VOLT_12','VOLT_5', 'VOLT_MAIN', 'REF_AVG','REF_STD','SW_DARK',...
                'SPEC_AVG','UV','CHECK'};
imosVarName =  {'NTRI','NTRI', 'AUX1','AUX2','AUX3','RMSe','TEMP_INT','TEMP_SPEC',...
                'TEMP_LAMP','LAMP_TIME','HUMIDITY','VOLT_12','VOLT_5', 'VOLT_MAIN', 'REF_AVG','REF_STD','SW_DARK',...
                'SPEC_AVG','UV',''};
            
nameMap = containers.Map(satlanticVar,imosVarName);

% output struct
sample_data            = struct;
sample_data.meta       = struct;
sample_data.variables  = {};
sample_data.dimensions = {};
% 
sample_data.toolbox_input_file          = filenames{1};
% 
% % The instrument_model field will be overwritten from the file header
sample_data.meta.instrument_make        = 'Satlantic';
sample_data.meta.instrument_model       = 'ISUS';
sample_data.meta.instrument_serial_no   = '';
sample_data.meta.featureType            = mode;

% calibration file data
i0_data = [];
eno3_data = [];
asw_data = [];
tasw_data = [];

% datafile scan lines
ts_data = cell(2,0);

cal_files = 1;
dat_files = 1;

% if we're given a .zip file, unzip to a temp directory
if (endsWith(filenames{1}, '.zip'))
    filenames = sort(unzip(filenames{1}, 'temp'));
end

% pre-process the xml files to get the frame structure
frameFormat = {};
fn_no = 1;
while (fn_no <= size(filenames,2))
    if (endsWith(filenames{fn_no}, '.xml'))
        fprintf('pre-processing : %s\n', filenames{fn_no});
    
        frameFormat = parseSatlanticXML(filenames{fn_no});
        sample_data.meta.instrument_model       = frameFormat.identifier;
        sample_data.meta.instrument_serial_no   = frameFormat.serialNumber;
    end
    fn_no = fn_no + 1;
end

% process the .CAL and .DAT files
fn_no = 1;
while (fn_no <= size(filenames,2))
    
    fprintf('processing : %s\n', filenames{fn_no});

    % read a .CAL file into a cal_instance
    if (endsWith(filenames{fn_no}, '.CAL'))

        [wlen, eno3, asw, tasw, i0] = readISUScalfile(filenames{fn_no});

        i0_data(cal_files,:) = i0;
        eno3_data(cal_files,:) = eno3;
        asw_data(cal_files,:) = asw;
        tasw_data(cal_files,:) = tasw;

        cal_files = cal_files + 1;

    end
    calInstance = 1:size(i0_data,1);

    % read data files into ts_data
    if (endsWith(filenames{fn_no}, '.DAT'))

        [datalines] = readISUSdatafile(filenames{fn_no});
        if (size(datalines,1) == 2) && (size(datalines,2) > 0)
            if (size(ts_data,2) > 0)
                if (size(ts_data{2},1) == size(datalines{2},1))
                    ts_data = horzcat(ts_data, datalines);
                end
            else
                ts_data = datalines;
            end
        end
        
        dat_files = dat_files + 1;

    end
    
    fn_no = fn_no + 1;
end

% convert ts_data into a vector array
data = cell2mat(ts_data(2,:));

% create static dimensions
sample_data.variables{end+1}.name           = 'TIMESERIES';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(1);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'LATITUDE';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'LONGITUDE';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'NOMINAL_DEPTH';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];

% if we have a cal file, add a CALIBRATION instance dimension
if cal_files > 1
    sample_data.dimensions{1}.name          = 'CALIBRATION_INSTANCE';
    sample_data.dimensions{1}.typeCastFunc  = str2func(netcdf3ToMatlabType('short'));
    sample_data.dimensions{1}.data          = sample_data.dimensions{1}.typeCastFunc(calInstance);
end

% if we have a data file, add a wavelength dimension for the UV variable
if (dat_files > 1)
    sample_data.dimensions{3}.name          = 'wavelength';
    sample_data.dimensions{3}.typeCastFunc  = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{3}.name, 'type')));
    sample_data.dimensions{3}.data          = sample_data.dimensions{3}.typeCastFunc(wlen);
end

fprintf('processed %d data files %d cal files %d times\n', dat_files-1, cal_files-1, size(ts_data, 2));

% dimensions definition must stay in this order : T, Z, Y, X, others;
% to be CF compliant

if (dat_files > 1)
    % order here so that plots plot these three first
    
    % TODO: how to translate the different frames into netCDF, where they
    % have different time vectors, for now they are just parsed with the
    % same format (asciiFrame full light) format
    for varIdx = 4:4
        varSelector = frameFormat.asciiFrame{varIdx}.varSelector;
        varNames = frameFormat.asciiFrame{varIdx}.varName;
    
        % create timestamps from Satlantic DATE (yeardoy), TIME (hour) data
        yeardoy = data(cell2mat(varSelector(strcmp(varNames, 'DATE'))),:);
        hour = data(cell2mat(varSelector(strcmp(varNames, 'TIME'))),:);
        year = fix(yeardoy/1000);
        day = rem(yeardoy,1000);
        timestamp = datenum(year, 1, 1) + (day - 1) + hour/24;

        sample_data.dimensions{2}.name          = 'TIME';
        sample_data.dimensions{2}.typeCastFunc  = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{2}.name, 'type')));
        sample_data.dimensions{2}.data          = sample_data.dimensions{2}.typeCastFunc(timestamp');

        % create the variables to hold the data
        for i=3:size(varNames)
            sName = char(varNames(i));
            if (sum(strcmp(sName, satlanticVar)) > 0)
            iName = nameMap(sName);
                if (length(iName) > 0)
                    if size(varSelector{i}, 2) == 1
                        sample_data.variables{end+1}.dimensions = [2];
                    else
                        sample_data.variables{end+1}.dimensions = [2 3];
                    end
                    sample_data.variables{end}.name         = iName;
                    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
                    sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(data(cell2mat(varSelector(i)),:)');
                end
            end
        end
    end
    
    % create a list of frame types, this is like the netCDF gridded index for the data    
    type_data = strcmp(ts_data(1,:), frameFormat.asciiFrame{1}.LineId); 
    for i=2:4
        type_data = type_data + strcmp(ts_data(1,:), frameFormat.asciiFrame{i}.LineId) * i; 
    end
    sample_data.variables{end+1}.dimensions = [2];
    sample_data.variables{end}.name         = 'SCAN_TYPE';
    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType('byte'));
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
    sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(type_data');
end


% create output variables to hold the calibration file(s) data
if (cal_files > 1)
    sample_data.variables{end+1}.dimensions = [1 3];
    sample_data.variables{end}.name         = 'I0';
    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';

    sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(i0_data);

    sample_data.variables{end+1}.dimensions = [1 3];
    sample_data.variables{end}.name         = 'ESW';
    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';

    sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(asw_data);

    sample_data.variables{end+1}.dimensions = [1 3];
    sample_data.variables{end}.name         = 'ETSW';
    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';

    sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(tasw_data);

    sample_data.variables{end+1}.dimensions = [1 3];
    sample_data.variables{end}.name         = 'ENO3';
    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';

    sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(eno3_data);

end
   
% plot calibration files
% figure(1); clf
% wlenIdx = getVar(sample_data.dimensions, 'wavelength');
% sh(1)=subplot(2,1,1);
% i0Idx = getVar(sample_data.variables, 'I0');
% plot(sample_data.dimensions{1,wlenIdx}.data, sample_data.variables{1,i0Idx}.data); grid on
% sh(2)=subplot(2,1,2); hold on
% i0Idx = getVar(sample_data.variables, 'ESW');
% plot(sample_data.dimensions{1,wlenIdx}.data, sample_data.variables{1,i0Idx}.data); grid on; 
% i0Idx = getVar(sample_data.variables, 'ENO3');
% plot(sample_data.dimensions{1,wlenIdx}.data, sample_data.variables{1,i0Idx}.data); grid on; 
% i0Idx = getVar(sample_data.variables, 'ETSW');
% plot(sample_data.dimensions{1,wlenIdx}.data, sample_data.variables{1,i0Idx}.data); grid on; ylim([0, 0.06]);
% legend('ESW', 'ENO3', 'ETSW');
% linkaxes(sh,'x');
% xlim([sample_data.dimensions{1,wlenIdx}.data(1) sample_data.dimensions{1,wlenIdx}.data(end)]);

end

% function to read ISUS calibration file, which holds the reference
% sprectra, and attenuation coefficents of sea-water, temperature*seawater
% and nitrate.
function [wlen, eno3, asw, tasw, i0] = readISUScalfile(filename)

    % Read a ISUS calibration file

    fid = fopen(filename); % Opens the calibration file

    cline = {}; % Creates an empty cell array to store the lines of data from the .dat file in

    i = 1;
    while ~feof(fid)
       tline = fgetl(fid);

        if strncmp(tline,'E',1) % only the extinction lines
            cline{i,1} = tline;
            comma = strfind(tline, ',');
            if (size(comma,2) == 5)
                ext(i,:) = sscanf(cline{i}(comma(1):end),',%f',6);
                i = i + 1;
            end
        end

    end

    wlen = ext(:,1);
    eno3 = ext(:,2); % ************** SWAPPED ENO3 and ASW, in KJ cal file
    asw = ext(:,3); 
    tasw = ext(:,4);
    i0 = ext(:,5);
    
    fclose(fid);

end

% RAW read the data into an array of vectors
function [data] = readISUSdatafile(filename)

    fid = fopen(filename); % Opens the data file, needs to be altered to loop through all .dat files in the directory

    % extract the data from the lines in the file
    
    lines = textscan(fid, '%s%[^\n]', 'Delimiter', ',');
    
    % mark all the header lines out
    %datalines = ~(strcmp(lines{1} ,'SATNHR') | (strcmp(lines{1} ,'SATAP6')));
    datalines = ~cellfun(@isempty, regexp(lines{1}(:), 'SAT[SN][LD]F\d\d\d\d'));

    datacount = sum(datalines);
    dl = lines{2}(datalines);
    dh = lines{1}(datalines);
    
    % datalines is a cell array with col 1 as the header (frame type), col 2 is the data
    data = cell(2, datacount);
    data(1,:) = dh;
    
    % TODO: check checksum of line
    dlCount = 1;
    for i=1:datacount
        v = textscan(char(dl(i)), '%f', 'Delimiter', ',');
        if (size(v{1},1) > 270) % is there a better way of doing this test?
            data(2,dlCount) = v;
            dlCount = dlCount + 1;
        end
    end
    data = data(:,1:(dlCount-1));

    fclose(fid);
    
end

% parse ISUS xml file
function frames = parseSatlanticXML(filename)

    xDoc = xmlread(filename);

    InstrumentPackage = xDoc.getElementsByTagName('Instrument');
    InstrumentPackageItem = InstrumentPackage.item(0);

    serialNumber = char(InstrumentPackageItem.getAttribute('serialNumber'));
    identifier = char(InstrumentPackageItem.getAttribute('identifier'));

    VarAsciiFrame = InstrumentPackageItem.getElementsByTagName('VarAsciiFrame');

    frames = {};

    frames.identifier = identifier;
    frames.serialNumber = serialNumber;

    for i=0:VarAsciiFrame.getLength()-1
        VarAsciiFrameItem = VarAsciiFrame.item(i);
        VarAsciiFrameItemIdentifier = VarAsciiFrameItem.getAttribute('identifier');
        % fprintf('identifier = %s\n', VarAsciiFrameItemIdentifier);
        frames.asciiFrame{i+1}.id = char(VarAsciiFrameItemIdentifier);

        SpectrumType = VarAsciiFrameItem.getElementsByTagName('SpectrumType');
        SpectrumTypeItem = SpectrumType.item(0);
        SpectrumTypeItemString = SpectrumTypeItem.getFirstChild().getData();

        Base = VarAsciiFrameItem.getElementsByTagName('Base');
        BaseItem = Base.item(0);
        BaseString = BaseItem.getFirstChild().getData();
        SerialNumber = VarAsciiFrameItem.getElementsByTagName('SerialNumber');
        SerialNumberItem = SerialNumber.item(0);
        SerialNumberString = SerialNumberItem.getFirstChild().getData();
        
        % fprintf('SpectrumType %s\n', SpectrumTypeItemString);
        frames.asciiFrame{i+1}.LineId = [char(BaseString) char(SerialNumberString)];
        frames.asciiFrame{i+1}.type = char(SpectrumTypeItemString);

        SensorFieldGroup = VarAsciiFrameItem.getElementsByTagName('SensorFieldGroup');
        for j=0:SensorFieldGroup.getLength()-1
            SensorFieldGroupItem = SensorFieldGroup.item(j);
            SensorFieldGroupName = SensorFieldGroupItem.getElementsByTagName('Name').item(0).getFirstChild().getData();
            % fprintf('%2d : %s\n', j, SensorFieldGroupName);
            frames.asciiFrame{i+1}.sensorGroup{j+1}.Name = char(SensorFieldGroupName);

            SensorField = SensorFieldGroupItem.getElementsByTagName('SensorField');
            for k=0:SensorField.getLength()-1
                SensorFieldItem = SensorField.item(k);
                SensorFieldIdentifier = SensorFieldItem.getElementsByTagName('Identifier').item(0).getFirstChild().getData();
                SensorFieldSequence = SensorFieldItem.getElementsByTagName('Sequence').item(0).getFirstChild().getData();
                % fprintf(' %2d : sequence %2s name : %s-%s\n', k, SensorFieldSequence, SensorFieldGroupName, SensorFieldIdentifier);            
                frames.asciiFrame{i+1}.sensorGroup{j+1}.sensorField{k+1}.Name = char(SensorFieldIdentifier);
                frames.asciiFrame{i+1}.sensorGroup{j+1}.sensorField{k+1}.Sequence = uint16(str2double((SensorFieldSequence)));
            end
        end
        
        % create a cell array with the names and data selector from the data lines
        frames.asciiFrame{i+1}.varName = cell(size(frames.asciiFrame{i+1}.sensorGroup,2)-2, 1);
        frames.asciiFrame{i+1}.varSelector = cell(size(frames.asciiFrame{i+1}.sensorGroup,2)-2, 1);
        for sensorGroup = 1:size(frames.asciiFrame{i+1}.sensorGroup,2)
            frames.asciiFrame{i+1}.varName{sensorGroup} = frames.asciiFrame{i+1}.sensorGroup{sensorGroup}.Name;
            selector = [];
            for sensorField = 1:size(frames.asciiFrame{i+1}.sensorGroup{sensorGroup}.sensorField,2)
                selector(sensorField) = frames.asciiFrame{i+1}.sensorGroup{sensorGroup}.sensorField{sensorField}.Sequence;
            end
            frames.asciiFrame{i+1}.varSelector{sensorGroup} = selector;
        end        
    end
end


