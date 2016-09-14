function [ outclass, f ] = svmclassifydist( svmStruct,sample, varargin )
%[outclass, f] = svmclassifydist( svmStruct, sample )


% set defaults
plotflag = false;

% check inputs
narginchk(2, Inf);

% deal with struct input case
if ~isstruct(svmStruct)
    error(message('stats:svmclassify:TwoInputsNoStruct'));
end

if ~isnumeric(sample) || ~ismatrix(sample)
    error(message('stats:svmclassify:BadSample'));
end

if size(sample,2)~=size(svmStruct.SupportVectors,2)
    error(message('stats:svmclassify:TestSizeMismatch'));
end

% deal with the various inputs
if nargin > 2
    if rem(nargin,2) == 1
        error(message('stats:svmclassify:IncorrectNumberOfArguments'));
    end
    okargs = {'showplot','-compilerhelper'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname)));
        if isempty(k)
            error(message('stats:svmclassify:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('stats:svmclassify:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % plotflag ('SHOWPLOT')
                    plotflag = opttf(pval,okargs{k}); 
                case 2 % help the compiler find required function handles by including svmtrain
                    svmtrain(eye(2),[1 0]);
            end
        end
    end
end

groupnames = svmStruct.GroupNames;

% check group is a vector -- though char input is special...
if ~isvector(groupnames) && ~ischar(groupnames)
    error(message('stats:svmclassify:GroupNotVector'));
end

% grp2idx sorts a numeric grouping var ascending, and a string grouping
% var by order of first occurrence
[~,groupString,glevels] = grp2idx(groupnames);  

% do the classification
if ~isempty(sample)
    % shift and scale the data if necessary:
    sampleOrig = sample;
    if ~isempty(svmStruct.ScaleData)
        for c = 1:size(sample, 2)
            sample(:,c) = svmStruct.ScaleData.scaleFactor(c) * ...
                (sample(:,c) +  svmStruct.ScaleData.shift(c));
        end
    end

    try
        [outclass,f] = svmdecision(sample,svmStruct);
    catch ME
        error(message('stats:svmclassify:ClassifyFailed', ME.message));
    end
    if plotflag

        if isempty(svmStruct.FigureHandles)
            warning(message('stats:svmclassify:NoTrainingFigure'));

        else
            try
                hAxis = svmStruct.FigureHandles{1};
                hLines = svmStruct.FigureHandles{2};
                hSV = svmStruct.FigureHandles{3};
                % unscale the data for plotting purposes
                [~,hClassLines] = svmplotdata(sampleOrig,outclass,hAxis); 
                trainingString = strcat(cellstr(groupString),' (training)');
                sampleString = strcat(cellstr(groupString),' (classified)');
                legend([hLines(1),hClassLines(1),hLines(2),hClassLines(2),hSV],...
                    {trainingString{1},sampleString{1},...
                    trainingString{2},sampleString{2},'Support Vectors'});
            catch ME
                warning(message('stats:svmclassify:DisplayFailed', ME.message));
            end
        end
    end
    outclass(outclass == -1) = 2;
    unClassified = isnan(outclass);
    outclass = glevels(outclass(~unClassified),:);
    if any(unClassified)
     
         try
             outclass = statinsertnan(unClassified,outclass);
         catch ME
             if ~isequal(ME.identifier,'stats:statinsertnan:LogicalInput')
                 rethrow(ME);
             else
                 error(message('stats:svmclassify:logicalwithNaN'));
             end
         end
    end

else
    outclass = [];
end


