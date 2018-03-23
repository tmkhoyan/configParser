function master()
% Masterfile to use config file

%   Author: Tigran Mkhoyan
%   Delft University of Technology, 2017




%% QO: Settings


configpath = 'config/config.txt';
%read config
[options] = readConfig(configpath,loadOptionsList,'//','{}',false); % varargin: 1=commentdelim, 2=headerdelim 3: warnings




    function optionsList = loadOptionsList
        %pass the options List; %caracters must match the txt file
        optionsList = {...
            'path'
            'paths'
            'boolean'
            'vector'
            'cell'
            };
    end

    
end