function sample_data = SBEcnvParse( fn, mode )
%SBEddParse Parse a raw '.cnv' file containing Sea Bird Data Conversion Output.
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

if iscellstr(fn)
    filename = fn{1};
else
    filename = fn;
end

if ~exist(filename, 'file')
    e = sprintf('File not found : %s\n', filename);
    error(e);
end

% save file size and open file; this will throw an error if file doesn't exist
filedir = dir(filename);
filesize = filedir.bytes;
fid = fopen(filename, 'rt');

T = readtable('Parser/seabirdNames.xlsx');

hardware_expr   = '* <HardwareData DeviceType=''(\S+)'' SerialNumber=''(\S+)''>';

name_expr       = '# name (\d+) = (.*): *(.*)';
nvalues_expr    = '# nvalues = (\d+)';
nquant_expr     = '# nquan = (\d+)';
end_expr        = '*END*';

use_expr        = '<(Use.*)>(.?)<\/\1>';
equa_group      = '<(.*) equation="(\d+)" >';

sensor_start    = '<sensor(.*)>';
sensor_end      = '</sensor(.*)>';
comment         = '<!--(.+?)-->';
tag             = '<(.+?)>(.+)</\1>';

header_expr  = '^\* Sea-Bird (.*?) *?Data File\:';

% output struct
sample_data            = struct;
sample_data.meta       = struct;
sample_data.variables  = {};
sample_data.dimensions = {};

% sample data and names, gleaned from the file
samples                = [];
varNames               = {};
    
sample_data.toolbox_input_file          = filename;

% The instrument_model field will be overwritten from the file header
sample_data.meta.instrument_make        = 'Sea-bird Electronics';
%sample_data.meta.instrument_model       = 'unknown';
%sample_data.meta.instrument_serial_no   = 'unknown';
sample_data.meta.featureType            = mode;

header = true; % start out looking for header
sensor = false;

lineno = 1;
datalineno = 1;
nSamples = 0;
nDataValues = 0;
sample_interval = 0;
dataNameFull = {};

% Read file 
nSensors = 0;
sensorN = 0;
sensorTag = 1;
sensorStruct = [];
use = -1;
equaNo = -1;
equaStr = '';

line = fgetl(fid);
while ischar(line)
    
    % parse the header
    if header
        if sensor
            % disp(line)
            % check for sensor expressions
            tkn = regexp(line, sensor_end, 'tokens');
            if ~isempty(tkn) 
                sensor = false;
                sensorN = 0;
            end
            % check for tags, save in sensor struct
            equa = regexp(line, equa_group, 'tokens');
            if ~isempty(equa) 
                equaNo = str2double(equa{1}{2});
                % fprintf("equation %d\n", equaNo);
                equaStr = equa{1}{1};
            end
            if regexp(line, ['</' equaStr '>'])
                equaNo = -1;
                eqaStr = '';
            end

            tkn = regexp(line, tag, 'tokens');
            if ~isempty(tkn) 
                useequ = regexp(line, use_expr, 'tokens');
                if ~isempty(useequ)
                    use = str2num(tkn{1}{2});
                    % fprintf('use %d\n', use);
                end
                if use == -1 || use == equaNo
                    % fprintf('tag %s = %s\n', tkn{1}{1}, tkn{1}{2});
                    if isempty(regexp(tkn{1}{1}, 'SerialNumber|CalibrationDate'))
                        tkn{1}{2} = str2double(tkn{1}{2}); % convert all others to double
                    end
                    sensorStruct(sensorN).tag(sensorTag) = tkn;
                    sensorTag = sensorTag + 1;
                end
            end
            % check for comment
            tkn = regexp(line, comment, 'tokens');
            if ~isempty(tkn) 
                % fprintf('comment %s\n', tkn{1}{1});
                if sensorN == 0 % only the first sensor
                    sensorN = nSensors;
                    sensorStruct(sensorN).Name = strtrim(tkn{1}{1});
                    use = -1;
                    equaNo = -1;
                end
            end
        else
            % check for sensor expressions
            tkn = regexp(line, sensor_start, 'tokens');
            if ~isempty(tkn) 
                sensor = true;
                nSensors = nSensors + 1;
                sensorTag = 1;
            end
            % check for name expressions
            tkn = regexp(line, name_expr, 'tokens');
            if ~isempty(tkn) 
                % create a vector of names and imos names for the variables
                dataValue = str2double(tkn{1}{1});
                idx = find(strcmp(tkn{1}{2}, T.ShortName), 1);
                if ~isempty(idx)
                    code(dataValue+1) = idx;         
                end
                dataNameFull(dataValue+1) = tkn{1}(3);
            end
            % check for number of values expressions
            tkn = regexp(line, nvalues_expr, 'tokens');
            if ~isempty(tkn) 
                nSamples = str2double(tkn{1}{1});            
                time = zeros(nSamples,1);
                samples = zeros(nSamples, nDataValues);
            end
            % check for number of quanties expressions
            tkn = regexp(line, nquant_expr, 'tokens');
            if ~isempty(tkn) 
                nDataValues = str2double(tkn{1}{1});            
            end
            % check for instrument type
            tkn = regexp(line, header_expr, 'tokens');
            if ~isempty(tkn) 
                sample_data.meta.instrument_model       = tkn{1};
            end
            % check for hardware
            tkn = regexp(line, hardware_expr, 'tokens');
            if ~isempty(tkn) 
                sample_data.meta.instrument_model       = tkn{1}{1};
                sample_data.meta.instrument_serial_no   = tkn{1}{2};
            end
            % every thing after the *END* is data
            tkn = regexp(line, end_expr, 'tokens');
            if ~isempty(tkn) 
                header = false;   
                samples = fscanf(fid, '%f', [nDataValues Inf])'; % readd the datablock after the END
        end
        end
    elseif ~header
        % read data line, using a , as separator
        %samples(datalineno,:) = cell2mat(textscan(line, '%f', nDataValues, 'Delimiter', '\t'));
        %samples = fscanf(fid, '%f', [nDataValues Inf]);
        %datalineno = datalineno + 1;
    end
    line = fgetl(fid);
    lineno = lineno + 1;
end

datalineno = size(samples, 1);
fprintf('samples read %d\n', datalineno);

% trim sample data to number of lines read
time = time(1:(datalineno-1));
samples = samples(1:(datalineno-1),:);

% massage any variables here as needed
tidx = find(strcmp(T.ShortName(code), 'timeK')); % elapsed time (seconds since 01-Jan-2000)
if ~isempty(tidx)
    time(:) = samples(:,tidx) / 86400 + datenum(2000,1,1);
end
tidx = find(strcmp(T.ShortName(code), 'timeY')); % elapsed time (seconds since 1970-01-01)
if ~isempty(tidx)
    time(:) = samples(:,tidx) / 86400 + datenum(1970,1,1);
end

% set the sample interval if we did not read it from the header
if (sample_interval < 0)
    sample_data.meta.instrument_sample_interval = median(diff(time*24*3600));
else
    sample_data.meta.instrument_sample_interval = sample_interval;
end

% copy the data into the sample_data struct
sample_data.dimensions{1}.name          = 'TIME';
sample_data.dimensions{1}.typeCastFunc  = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{1}.name, 'type')));
sample_data.dimensions{1}.data          = sample_data.dimensions{1}.typeCastFunc(time);

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

%fprintf('Var Names:\n');
%disp(T.FullName(code(:)));
%disp(dataNameFull');

sensorNames = {};
for i = 1:size(sensorStruct,2)
    sn = regexp(sensorStruct(i).Name, '[^,]*, (.*)', 'tokens');
    if ~isempty(sn)
        sensorNames{i} = sn{1}{1};
        sensorVars(i,:) = startsWith(T.FullName(code(:)), sensorNames(i));
    end
end

% Create sample data for each variable in the file
for k = 1:length(code)

  varName = T.IMOSCode(code(k));
  % variables not to put into Variables section
  if strcmp(varName, 'TIME')
      continue;
  end

  % create a new variable
  sample_data.variables{end+1}.dimensions = 1;

  % get the sea-bird table units
  unit = T.Units(code(k));

  % if there is no IMOS code, use the sea-bird Short Name from the table
  if isempty(varName{1})
      varName = T.ShortName(code(k));
      varName = regexprep(varName, '[/-]', '_');
      varName = strrep(varName, 'é', 'u');
      sample_data.variables{end}.unit = unit{1};
  else
      sample_data.variables{end}.unit = imosParameters(varName{1}, 'unit');
  end
      
  % dimensions definition must stay in this order : T, Z, Y, X, others;
  % to be CF compliant
  sample_data.variables{end}.name         = varName{1};
  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType('float')); % imosParameters(varName{1}, 'type')
  sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';

  sample_data.variables{end}.sea_bird_unit         = unit{1};

  %fullName = T.FullName(code(k));
  %sample_data.variables{end}.sea_bird_full_name         = fullName{1};
  sample_data.variables{end}.file_full_name         = dataNameFull{k};
  
  % add any tags from the processing for the sensors for the variables
  sensors = sensorVars(:,k);
  if sum(sensors) ~= 0
      sensor_tags = sensorStruct(sensors).tag;

      for sn = 1:size(sensor_tags,2)
          cal = ['calibration_' sensor_tags{sn}{1}];
          sample_data.variables{end}.(cal)         = sensor_tags{sn}{2};
      end
  end
  sample_data.variables{end}.data = sample_data.variables{end}.typeCastFunc(samples(:,k));
  
  if (strcmp(varName, 'PRES_REL'))
        % let's document the constant pressure atmosphere offset previously
        % applied by SeaBird software on the absolute presure measurement
        sample_data.variables{end}.applied_offset = sample_data.variables{end}.typeCastFunc(-14.7*0.689476);
  end
end
