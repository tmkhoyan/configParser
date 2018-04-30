
function [configOptionsStruct, configOptionsCell, allheaders] = readConfig(filename,varargin)
% Config parser for ASCII based files

%   Author: Tigran Mkhoyan
%   Delft University of Technology, 2017

% see below for description
headerDelim = [];
maxheaders = [];
structnamefieldfillelemn = [];
setOptargs;

fileID = fopen(filename);
charsizecomment = numel(commentDelim); %#ok<SHVAI>

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
        
        if (nextline([1 end]) == headerDelim)  % % new header
            
            rowcount = 2;
            configOptionsCell(1,headercount) = {nextline};    %header
            allheaders(headercount,1) = {nextline(2:end-1)}; %row two is header without the brakets
            headercount = headercount + 1;
            % optionList
        elseif (headercount>1)
            if (isempty(str2num(nextline))) %#ok<*ST2NM> %%&& ~strcmp(nextline,'[]'))
                %if (isempty(str2double(nextline))) %%&& ~strcmp(nextline,'[]'))
                configOptionsCell(rowcount,headercount-1) = {nextline}; %option
                noncharIdx(rowcount,headercount-1) = 0;
            else
                configOptionsCell(rowcount,headercount-1) = {str2num(nextline)}; %convert to numbers
                %configOptionsCell(rowcount,headercount-1) = {str2double(nextline)}; %convert to numbers
                noncharIdx(rowcount,headercount-1) = 1;
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
linenum = loopcount; %#ok<NASGU>
emptyRows = ~sum(~cellfun('isempty',configOptionsCell),2); %emty rows
configOptionsCell(emptyRows,:) = [];
noncharIdx(emptyRows,:) = [];
% checkheaders = configOptionsCell{~noncharIdx};
% checkvals    = configOptionsCell{2,k};

noncharIdx(1,:) = 1; % dont covert headers
noncharIdx = logical(noncharIdx);
allheaders = allheaders(1:size(configOptionsCell,2),:);

configOptionsCell(cellfun('isempty',configOptionsCell)) = {'~NaN'}; % fill all empty fields with nan

%TODO consider using function handle instead of parsing pwd. However then
%option must be accesed as opt.myopt();
%mydir = @() pwd;
%configOptionsCell(strcmp(configOptionsCell,'pwd')) = {mydir}; //to make
%it an aninymous function. However in usage option must be called with
%opt.mydir()
configOptionsCell(strcmp(configOptionsCell,'pwd')) = {pwd()};

%Update2.0: functionality header addition:
%c1 = configOptionsCell;

for k = 1:numel(configOptionsCell(1,:))
    checkheader = configOptionsCell{1,k};
    checkval    = configOptionsCell{2,k}; %TODO:takes only first value if the header has multiple values wont work
    if(ischar(checkval))
        configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),checkheader,checkval);
    end
end
%check for double path //

configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),'//','/');

% %in case [] is config option make sure its numeric not string
configOptionsCell(strcmp(configOptionsCell,'[]')) = {[]};
%create struct

%Update 2.2: recognize = sign and split header options
% for k = 1:numel(configOptionsCell(1,:))
%     %for m =1:numel(configOptionsCell())
%     checkval    = configOptionsCell{2,k}; %TODO:takes only first value if the header has multiple values wont work
%     if(ischar(checkval))
%         configOptionsCell(~noncharIdx) = regexprep(configOptionsCell(~noncharIdx),'=','split');
%     end
% end

%generate temp option cell and replace non char values
tmpconfigOptionsCell = configOptionsCell;
tmpconfigOptionsCell(noncharIdx) = {''};

dynamicHeaderValIdx = find(~~cellfun(@sum,(regexp(tmpconfigOptionsCell,'=')))); % get the ids of cells to be replaced

toSplitIdx = ~~(zeros(size(noncharIdx)));
toSplitIdx(dynamicHeaderValIdx) = 1;

configOptionsCell(toSplitIdx) = strtrim(regexp(configOptionsCell(toSplitIdx),'=','split'));
toSplitHeaderIdx = sum(toSplitIdx)>=1; %indeces of headers with special = operator

%zz = regexp(configOptionsCell(~noncharIdx),'=','split');
%zz2 = regexp(configOptionsCell(:,9),'=','split');
% update: check if header contains same variable twice remove (only first one is used)
% [~,uniqueindeces] = unique(allheaders(:,1));
% uniqueindeces = sort(uniqueindeces); % keep order
% allheaders = allheaders(uniqueindeces,:);
% headercount = numel(uniqueindeces)+1;
% configOptionsCell = configOptionsCell(:,uniqueindeces);

%UPDATE for simulink parser
%replace headers names with valid fielnames for structure
charToreplace = {'/',' '};
configOptionsCell(1,:) = regexprep(configOptionsCell(1,:),charToreplace,structnamefieldfillelemn);
allheaders(:,1) = regexprep(allheaders(:,1),charToreplace,structnamefieldfillelemn);

[allheaders, uniqueindeces] = checkUniqness(allheaders);
headercount = numel(uniqueindeces)+1;
configOptionsCell = configOptionsCell(:,uniqueindeces);
%do the same for
%match options if given


configOptionsStruct = struct();
if isempty(optionsList)
    optionsList = allheaders(:,1); %use all headers in case no specific headers are provded to build the struct
else
    %replace headers names with valid fielnames for structure
    optionsList = checkUniqness(optionsList);
    configOptionsCell(1,:) = regexprep(configOptionsCell(1,:),charToreplace,structnamefieldfillelemn);
end

%  configOptionsStruct = initStruct(optionsList);
optionsSelected = cell(max([allheaders{:,2}]+1),numel(optionsList));

if (headercount-1)>numel(optionsList)
    if (warnEnabled) %#ok<SHVAI>
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
    % fieldname(isspace(fieldname)) = structnamefieldfillelemn; already
    % done by regexprep
    if (max(rowsmatched)>2 || toSplitHeaderIdx(i)) % if we have special headers '=' put them al in cell container
        configOptionsStruct = setfield(configOptionsStruct,fieldname,configOptionsCell(rowsmatched(rowsmatched>=2),columnmatched)); %#ok<SFLD> %fieldvalue is cell
    else
        configOptionsStruct = setfield(configOptionsStruct,fieldname,configOptionsCell{rowsmatched(rowsmatched>=2),columnmatched}); %#ok<SFLD> %fieldvalue is just value
    end
    
end

configOptionsCell = optionsSelected;


%     function s = initStruct(fieldnames)
%     zz = fieldnames;
%     zz(:,2) = {[]};
%     zz = zz';
%     s = struct(zz{:});
%
%     end

    function [incell, uniqueIdx] = checkUniqness(incell)
        [~,uniqueIdx] = unique(incell(:,1)); % checks only first column
        uniqueIdx = sort(uniqueIdx); % keep order
        incell = incell(uniqueIdx,:); % returns full cell
    end

%UPDATE simulinkparser: now the runtimesettings can be set using two
%methods: 1. via providing ordered variable arguments. 2. Via providing
%name tag of runtime option followed by the value e.q.
%readconfig(...,'structnamefieldfillelemn','__',...) irrespective of the order of variable argument this usage of readConfig will pick the right configuration for the explicit setting 
    function setOptargs
        numvarargs  = length(varargin);
        
        optargs = {cell(0), '//', '{}',true, '_', 100};
        %optargs{1:numvarargs} = varargin;
        %UPDATE simulink parser. Now one option is parced
        runtimeSettings = {'optionsList', 'commentDelim', 'headerDelim', 'warnEnabled', 'structnamefieldfillelemn', 'maxheaders'};
        idx = cell(size(runtimeSettings));
        idxval = idx;
        id = idx;
        %idxval = ~~(zeros(numvarargs,numel(runtimeSettings)));
        
        tmp = (zeros(size(runtimeSettings)));
        % check if explicit options are provided
        for n=1:numel(runtimeSettings)
            tmp(n) = 1;
            idx{n} = strcmpi(varargin,runtimeSettings{n})';
            idxval{n} = circshift(idx{n},1);
            id{n} = (tmp*sum(idx{n})==1)'; % (tmp*(sum(idx{n})==1))'; gives doubles array!
            tmp = tmp*0;
         end
        
        if sum(cellfun(@sum,idx)>1) %if duplicates
            str = strcat(runtimeSettings,',');
            error('runtime setting ''%s'' specified more than once',[str{sum(zz,1)>1}])
            
        elseif ~(sum(cellfun(@sum,idx))) % regular way of setting runtimesettings
            % set defaults for optional inputs
            if numvarargs > 5
                error('functions:randRange:TooManyInputs', ...
                    'requires atmost 2 optional input');
            end
            %TODO: make this work for all variable arguments
            [optargs{1:numvarargs}] = varargin{:};
            [optionsList, commentDelim, headerDelim, warnEnabled, structnamefieldfillelemn, maxheaders] = optargs{:};
        else % parce from strings order does not matter
            %set optional argumenst and apliy specific runtime settings for the remaining
            %idxval = circshift(idx,1); %value must be the enry afer
            for n=1:numel(runtimeSettings)
                [optargs{id{n}}] = varargin{idxval{n}};
            end
            [optionsList, commentDelim, headerDelim, warnEnabled, structnamefieldfillelemn, maxheaders] = optargs{:};
            
        end
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