function [configOptionsStruct, configOptionsCell, allheaders] = readConfig(filename,varargin)
% Config parser for ASCII based files 

%   Author: Tigran Mkhoyan
%   Delft University of Technology, 2017

% see below for description

setOptargs;

fileID = fopen(filename);
charsizecomment = numel(commentDelim);

fgets(fileID); %skip header
N = 1000; %lines

configOptionsCell = cell(N/10,2);
allheaders = cell(N/10,2);

headercount = 1;
loopcount = 1;
rowcount = 2; %to make sure first heder gets

noncharIdx = zeros(size(configOptionsCell));
%put back old
while (~feof(fileID) && loopcount < N && ((headercount+rowcount) <= (maxheaders+10)))
    
    %read
    nextline = strtrim(fgetl(fileID));
    numchar = numel(nextline);
    if (~isempty(nextline) && sum(nextline(1:min(numchar,charsizecomment)) == commentDelim) ~= charsizecomment)
        
        if (nextline([1 end]) == headerDelim)  % new header
            
            rowcount = 2;
            configOptionsCell(1,headercount) = {nextline};    %header
            allheaders(headercount,1) = {nextline(2:end-1)}; %row two is header without the brakets
            headercount = headercount + 1;
            % optionList
        elseif (headercount>1)
            if (isempty(str2num(nextline))) %%&& ~strcmp(nextline,'[]'))
            configOptionsCell(rowcount,headercount-1) = {nextline}; %option
            else 
            configOptionsCell(rowcount,headercount-1) = {str2num(nextline)}; %convert to numbers    
            noncharIdx(rowcount,headercount-1) = [1];
            end 
            rowcount = rowcount + 1;
            allheaders(headercount-1,2) = {rowcount-2}; %update number of elements
        end
    end
    
    loopcount = loopcount + 1;
    
    %TODO fix the number of headers part
    if (headercount+rowcount) == (maxheaders+10)
    warning('Maxixmum number of header elements reached. Ignoring extra elemenst. Use read config with optional argumens to increase')
    end
end

%cleanup options
linenum = loopcount;
emptyRows = ~sum(~cellfun('isempty',configOptionsCell),2); %emty rows
configOptionsCell(emptyRows,:) = [];
noncharIdx(emptyRows,:) = [];
noncharIdx(1,:) = 1; % dont covert headers
noncharIdx = logical(noncharIdx);
allheaders = allheaders(1:size(configOptionsCell,2),:);

configOptionsCell(cellfun('isempty',configOptionsCell)) = {'~NaN'}; % fill all empty fields with nan

 mydir = @() pwd;
 %configOptionsCell(strcmp(configOptionsCell,'pwd')) = {mydir}; //to make
 %it an aninymous function. However in usage option must be called with
 %opt.mydir()
 configOptionsCell(strcmp(configOptionsCell,'pwd')) = {pwd()};
 
 % adding variables together
 %c1 = configOptionsCell;
 for k = 1:numel(configOptionsCell(1,:))
     checkheader = configOptionsCell{1,k};
     checkval    = configOptionsCell{2,k}; %takes only first value if the header has multiple values wont work
     if(ischar(checkval))
        configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),checkheader,checkval);
     end
 end
 
 %check for double path //

  configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),'//','/');

 
 
% %in case [] is config option make sure its numeric not string
 configOptionsCell(strcmp(configOptionsCell,'[]')) = {[]};
%create struct 




%match options if given
configOptionsStruct = struct();
if isempty(optionsList)
    optionsList = allheaders(:,1); %use all headers in case no specific headers are provded to build the struct
end 
   
  %  configOptionsStruct = initStruct(optionsList);
    optionsSelected = cell(max([allheaders{:,2}]+1),numel(optionsList));
    
    if (headercount-1)>numel(optionsList)
        if (warnEnabled)
        warning(['Extra heading entries in config file, ',filename, ', ignored'])
        end
    end
    
    for i=1:numel(optionsList)
        
        columnmatched = strcmpi(allheaders(:,1),optionsList(i));
        if ~columnmatched
            error(['config option, ', optionsList{i}, ', not found'])
        else
            rowsmatched = 1:(1+allheaders{columnmatched,2});
        end
        optionsSelected(rowsmatched,i) = configOptionsCell(rowsmatched,columnmatched);
        
        %set struct fields
        fieldname  = optionsList{i};
        fieldname(isspace(fieldname)) = structnamefieldfillelemn;
        if max(rowsmatched)>2
        configOptionsStruct = setfield(configOptionsStruct,fieldname,configOptionsCell(rowsmatched(rowsmatched>=2),columnmatched)); %fieldvalue is cell
        else 
        configOptionsStruct = setfield(configOptionsStruct,fieldname,configOptionsCell{rowsmatched(rowsmatched>=2),columnmatched}); %fieldvalue is just value
        end 
        
    end


configOptionsCell = optionsSelected;


    function s = initStruct(fieldnames)
    zz = fieldnames;
    zz(:,2) = {[]};
    zz = zz';
    s = struct(zz{:});
    
    end

    function setOptargs
        numvarargs  = length(varargin);
        
        % set defaults for optional inputs
        if numvarargs > 5
            error('functions:randRange:TooManyInputs', ...
                'requires atmost 2 optional input');
        end
        
        optargs = {cell(0), '//', '{}',true, '_', 100};
        %optargs{1:numvarargs} = varargin;
        [optargs{1:numvarargs}] = varargin{:};
        [optionsList, commentDelim, headerDelim, warnEnabled, structnamefieldfillelemn, maxheaders] = optargs{:};
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
       
% Usage of 2.0 functionality: 
% 
% - as of 2.0 the parser will recognize pwd command. so using: 
% {mydir} --> {mydir}
% pwd     --> path/to/location/calling/master/function
% In other words if you call readConfig with mast.m located at ~/somedirectory/ the options.mydir = ~/somedirectory/
% 
% Parser will now also allow variable additions. In other words if you define a variable {mydir} you can use this to construct a new variable composed partially or fully of {mydir}. Example:
% 
% define:
%  {mydir}
% ~/somepath/
% use {mydir} to definecomposite path variable 
% 
% {myotherdir}
% {mydir}/../ 
% 
% Now the {myotherdir} will point to a folder level higher than {mydir} (this is due to addition ../)
% This can also be use for addition of any string variables. Note that the functionality doesn't work yet for numeric values. In case there is interest i might consider adding that.  