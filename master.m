function master()

% Masterfile to use config file  for ASCII based files 

%   Author:         Tigran Mkhoyan (email: t.  mkhoyan@tudelft.nl)
%                   Delft University of Technology, 2017


% ---------------------------------------------------------------------------------:

% Usage

configpath = 'config/config.txt';
%read config
[options, optionsCell] = readConfig(configpath,loadOptionsList,'//','{}',false); % varargin: 1=commentdelim, 2=headerdelim 3: warnings


%display 
optionsList = loadOptionsList;
T = cell2table(optionsCell');
varnames = circshift(T.Properties.VariableNames,1); varnames{1} = 'Headers';
T.Properties.VariableNames = varnames;  
%T = table(optionsList,optionsList,'Variablenames',{'variables' 'value'});
disp('   ----------------config options-----------------')
disp(' ')
disp(T)

% this is a local function to create the header names. 

    function optionsList = loadOptionsList
        %pass the options List; %caracters must match the txt file
        optionsList = {...
            'path'
            'paths'
            'boolean'
            'vector'
            'cell'
            'currentdir'
            'projectdir'
            'appendeddir'
            'correctedpath'
            };
    end

    
end

% ---------------------------------------------------------------------------------:

% Description/Usage:
% -A simle config parser for types of string,boolean, vector,cell and matrix values. 
%   The parser uses simple headers with customizable header and comment delimiters. 
%   In other words use a delimiter of your preference for comments ( e.g. // or %, or  ) and setting headers (e.g. {myoption1} or [myoption2]). 
%   De delimiters can be provided to readConfig as second and third input for comments and headers respectively. First option is the path to the config file.
% -The parser will generate a stuct with settings that are accecible via same header names. Example: 

%   in config file  --> set option with 
%   {myoption1}
%   1 2 3
%   in matlab file  --> acces vector [1 2 3] opt struct with
%   opt.myoption1
  
%   The parser will automatically recognize the following types: 
  
%  -(single) path/string    = string with pathname/caracter
%  -(multiple) path/        = 1xm cell with m pathnames listed after {header}
%  - boolean                = false or true logical (case sensitive!)
%  - vector                 = simple 1xm vector with m space separated values after header
%  - cell/matrix            = 2D matrix with m columns (space separated values) and n rows (lines after header). Note that the           function returns a cell. To convert to array simple do cell2mat(opt.mycellarray)
 
  
% -File type is also customizable (i.e. txt or conf or whatever)

% Usage: 
% -create a txt file similar to config.txt
%   e.g.
% -put desired header names in a cell with strings (order is not relevant but caracters must match) 
% -acces config with:
%     [options] = readConfig('path/to/config',headernames,'//','{}',false); % varargin: 1=commentdelim, 2=headerdelim 3: -
% -use additional options to customyze to needs
% loadOptionsList is a local function to create the header names. If this is not desired simply use:
       % optionsList = {...
       %      'path'
       %      'paths'
       %      'boolean'
       %      'vector'
       %      'cell'
       %      }; 