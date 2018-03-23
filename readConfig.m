function [configOptionsStruct, configOptionsCell, allheaders] = readConfig(filename,varargin)
% Config parser for ASCII based files 

%   Author: Tigran Mkhoyan
%   Delft University of Technology, 2017
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
allheaders = allheaders(1:size(configOptionsCell,2),:);

configOptionsCell(cellfun('isempty',configOptionsCell)) = {'~NaN'}; % fill all empty fields with nan

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

