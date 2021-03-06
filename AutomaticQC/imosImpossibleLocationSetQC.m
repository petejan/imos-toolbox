function [sample_data, varChecked, paramsLog] = imosImpossibleLocationSetQC( sample_data, auto )
%IMOSIMPOSSIBLELOCATIONSETQC Flags impossible Latitude and Longitude values 
% using IMOS sites information from imosSites.txt.
%
% Impossible location test described in Morello et al. 2011 paper.
%
% Inputs:
%   sample_data - struct containing the entire data set and dimension data.
%   auto        - logical, run QC in batch mode
%
% Outputs:
%   sample_data - same as input, with QC flags added for variable/dimension
%                 data.
%   varChecked  - cell array of variables' name which have been checked
%   paramsLog   - string containing details about params' procedure to include in QC log
%
% Author:       Guillaume Galibert <guillaume.galibert@utas.edu.au>
%

%
% Copyright (C) 2017, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.
% If not, see <https://www.gnu.org/licenses/gpl-3.0.en.html>.
%

narginchk(1, 2);
if ~isstruct(sample_data), error('sample_data must be a struct'); end

% auto logical in input to enable running under batch processing
if nargin<2, auto=false; end

varChecked = {};
paramsLog  = [];

qcSet    = str2double(readProperty('toolbox.qc_set'));
passFlag = imosQCFlag('good',           qcSet, 'flag');
failFlag = imosQCFlag('probablyBad',    qcSet, 'flag');

idLon = getVar(sample_data.variables, 'LONGITUDE');
dataLon = sample_data.variables{idLon}.data;
idLat = getVar(sample_data.variables, 'LATITUDE');
dataLat = sample_data.variables{idLat}.data;

if ~isempty(dataLon) && ~isempty(dataLat)
    % get details from this site
%     site = sample_data.meta.site_name; % source = ddb
%     if strcmpi(site, 'UNKNOWN'), site = sample_data.site_code; end % source = global_attributes file
    site = sample_data.site_code;
    
    site = imosSites(site);
    
    % test if site information exists
    if isempty(site)
        if isempty(sample_data.site_code)
            fprintf('%s\n', ['Warning : ' 'No site code documented to '...
                'perform impossible location QC test']);
        else
            fprintf('%s\n', ['Warning : ' 'No site information found for ' sample_data.site_code ' in imosSites.txt to '...
                'perform impossible location QC test']);
        end
    else
        lenData = length(dataLon);
        
        % initially all data is bad
        flagLon = ones(lenData, 1, 'int8')*failFlag;
        flagLat = flagLon;
        
        %test location
        if isnan(site.distanceKmPlusMinusThreshold)
            paramsLog = ['longitudePlusMinusThreshold=' num2str(site.longitudePlusMinusThreshold) ...
                ', latitudePlusMinusThreshold=' num2str(site.latitudePlusMinusThreshold)];
            
            % test each independent coordinate on a rectangular area
            iGoodLon = dataLon >= site.longitude - site.longitudePlusMinusThreshold && ...
                dataLon <= site.longitude + site.longitudePlusMinusThreshold;
            iGoodLat = dataLat >= site.latitude - site.latitudePlusMinusThreshold && ...
                dataLat <= site.latitude + site.latitudePlusMinusThreshold;
        else
            paramsLog = ['distanceKmPlusMinusThreshold=' num2str(site.distanceKmPlusMinusThreshold)];
            
            % test each couple of coordinate on a circle area
            if (site.latitude == dataLat) && (site.longitude == dataLon)
                obsDist = 0; % case the actual position is exactly the nominal location!
            else
                obsDist = WGS84dist(site.latitude, site.longitude, dataLat, dataLon);
            end
            
            iGoodLon = obsDist/1000 <= site.distanceKmPlusMinusThreshold;
            iGoodLat = iGoodLon;
        end
        
        if any(iGoodLon)
            flagLon(iGoodLon) = passFlag;
        end
        
        if any(iGoodLat)
            flagLat(iGoodLat) = passFlag;
        end
        
        sample_data.variables{idLon}.flags = flagLon;
        sample_data.variables{idLat}.flags = flagLat;
        
        varChecked = {'LATITUDE', 'LONGITUDE'};
        
        if any(~iGoodLon) || any(~iGoodLat)
            error('Impossible location QC test failed => Check deployment database values for latitude/longitude or thresholds in imosSites.txt!');
        end
    end
end