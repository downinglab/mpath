clear
clc
close all
% January 16, 2024: Principal Component Analysis for Long-Read WGBS Data
% Goal: Develop robust method of analyzing multiple metrics (Uniformity,
% Pearson Corr, etc.) and their ability to separate between Nascent/Mature
% DNA reads

% Input: .txt files for Nascent and Mature datasets; Columns correspond to
% Read Length, WGBS fraction, and generated metrics; rows each represent
% different reads.

%ensure that the two input tables were generated from the same script
%(Python or MATLAB versions)

%% file name retrieval
[~,~,filen_Nasc] = browseFile("nascent");
[~,~,filen_Mature] = browseFile("mature");

%% reading tables
TT=readtable(filen_Nasc,"VariableNamingRule","preserve");
TT2=readtable(filen_Mature,"VariableNamingRule","preserve");

%% reformat tables if needed
TT = reformatTable(TT);
TT2 = reformatTable(TT2);
Zmark=size(TT,1); % This mark checks how many reads are in the Nascent datafile

%
%% unwanted cell arrays of character arrays to double
%the following two lines to set metric columns back to
%double precision arrays from character arrays
%the cellCheck function assumes that the first column in the table is the
%read ID and the rest should be numerical data
TT = [TT(:,1), varfun(@(r) cellCheck(r),TT(:,2:end))];
TT2 = [TT2(:,1), varfun(@(r) cellCheck(r),TT2(:,2:end))];

%% Concatenate Nascent and Mature read data
%track which columns are not named with the second arg
getDataColumnIndicesByColumnNames = ~contains(TT.Properties.VariableNames,["ReadID","ReadLength","read_id"]);
%track which columns are of class double
getDataColumnIndicesByType = table2array(varfun(@(r) matches(class(r),'double'),TT));
%get the intersection of the two conditions
getDataColumnIndices = getDataColumnIndicesByColumnNames & getDataColumnIndicesByType;
%if error "tabular/vertcat" appears, ensure that the nascent and mature
%tables are generated from the same script.
TTT=table2array([TT(:,getDataColumnIndices);TT2(:,getDataColumnIndices)]);
%% pca analysis
% Perform PC Analysis (Simple)
[coeff,score,~,~,explained,mu]=pca(TTT);

% Meta Params
fonttt=24;

%% append nascent and mature tables with pca scores
nascentPCA1Scores = score(1:Zmark,1);
nascentPCA2Scores = score(1:Zmark,2);
maturePCA1Scores = score(Zmark+1:end,1);
maturePCA2Scores = score(Zmark+1:end,2);
TT = addvars(TT,nascentPCA1Scores,nascentPCA2Scores);
TT2 = addvars(TT2,maturePCA1Scores,maturePCA2Scores);

%% write the new tables to file with appended PCA1,2 scores
%comment the open block comment below if you'd like to output the input data
%tables with the PCA 1,2 scores appended and write as new tables.
%{
selectedPath = string(uigetdir([],"Select file save destination folder"));
[~,fileNameNasc,fileExtNasc] = fileparts(filen_Nasc);
[~,fileNameMature,fileExtMature] = fileparts(filen_Mature);
newNascFilePath = selectedPath + filesep + fileNameNasc + ...
    "_pca_scores_appended" + fileExtNasc;
newMatureFilePath = selectedPath + filesep + fileNameMature + ...
    "_pca_scores_appended" + fileExtMature;
writetable(TT,newNascFilePath);
writetable(TT2,newMatureFilePath);
%}

%% Visualization
figure(1)
imagesc(coeff)
%colorbar
title('Coefficient Matrix')
ax=gca;
ax.FontSize=fonttt;
ax.TickDir='out';
ax.TickLength=[0.05 0.05];
grid on
grid minor

figure(2)
plot(1:numel(explained),explained,'ok','LineWidth',1,'MarkerFaceColor','black')
title('Explained')
xlabel('PC No.')
ylabel('% of Variance')
yticks([0 20 40])
ax=gca;
ax.FontSize=fonttt;
grid minor

figure(3)
scatter(score(:,1),score(:,2),'filled','MarkerFaceColor','#4DBEEE','DisplayName','16hr')
hold on
scatter(score(1:Zmark,1),score(1:Zmark,2),'filled','MarkerFaceColor','[1 0.5 0.8]','DisplayName','0hr')
hold off
xlabel('PCA-1')
ylabel('PCA-2')
legend('Location','best')
ax=gca;
ax.FontSize=fonttt;
grid minor
grid on

%scatter histogram below

figure(4)
groupNameCell = cell(1,size(score,1));
groupBool = [true(1,Zmark),false(1,length(score(Zmark+1:end,1)))];
groupNameCell(groupBool) = {'0hr'};
groupNameCell(~groupBool) = {'16hr'};
scatterhistogram(score(:,1),score(:,2),'GroupData',groupNameCell,'HistogramDisplayStyle','smooth');
xlabel('PCA-1')
ylabel('PCA-2')
xlim([-1.75 1.75]);
ylim([-1.75 1.75]);

%% local functions
%convert a given table column expected to be data to type double if it is
%not already
function outputTableColumn = cellCheck(tableColumn)
    if (matches(class(tableColumn),'cell'))
        if (matches(class(tableColumn{1}),'char'))
            outputTableColumn = double(string(tableColumn));
        end
    elseif ~matches(class(tableColumn),'double')
        outputTableColumn = double(tableColumn);
    else
        outputTableColumn = tableColumn;
    end
end

%opens file dialog for user to select table text file
function [fileNames, folderPath, getDirectoryList] = browseFile(messageToDisplay)
    switch messageToDisplay
        case "nascent"
            messageToDisplay = "Select Nascent Data File";
        case "mature"
            messageToDisplay = "Select Mature Data File";
    end
    %both the file and path names are stored. extension and file names are
    %wildcards.
    [fileNames, folderPath] = uigetfile('*.*',messageToDisplay);
    %the full directory path is generated by concatenating the path and
    %file strings
    getDirectoryList = strcat(folderPath, fileNames);
end

%this function determines if the table was generated from the python or
%matlab versions of the metrics code. the two scripts generate tables with
%different formats where the python version will return an error when
%loaded into script "NascMatur_PCA_scatterhistogram.m". Therefore, table
%will be rearranged into a format recognized by the PCA script. 

%NOTE that the input tables need to be of the same dimension column-wise,
%or in other words, generated from the same script (python or MATLAB),
%otherwise an error will appear.
function [rearrangedTable] = reformatTable(currentTable)
    %get the current table's variable names
    getTableVariableNames = string(currentTable.Properties.VariableNames);
    %the expected table variable names from the python version of the
    %metrics script
    expectedPythonColumnNames = ["read_id","read_wgbs_distance","read_meth_ratio","distance_bin","num_pairs","smc","uniformity","pearson_r","pearson_p"];
    %if the variable names do not match that of the python generated table,
    %it can be assumed that the table came from the MATLAB metrics script.
    %currently not accounting for any other cases.
    if ~all(matches(getTableVariableNames,expectedPythonColumnNames))
        rearrangedTable = currentTable;
        return
    end
    %exclude the "distance_bin" and "num_pairs" columns since they will be
    %integrated into the updated variable names and not relevant to
    %analysis, respectively.
    getTableVariableNames = getTableVariableNames([1:3,6:9]);
    %expected bin names from the python metrics script
    binNames = ["bin_all", "bin_0", "bin_1", "bin_2", "bin_3", "bin_closest"];
    %get the table variable column that contains the bin identifier
    getBinColumn = string(table2array(currentTable(:,4)));
    %get the table variable column that contains the read IDs
    getReadID = string(table2array(currentTable(:,1)));
    %get the unique read IDs
    [uniqReadID, ~, ic] = unique(getReadID,"stable");
    %preallocate the new data table. 4 bin dependent data columns with 6
    %different bin conditions plus the two non bin dependent data columns
    %(read meth ratio and wgbs distance). change "4" in second argument to
    %5 if keeping the number of pairs column in the new table
    updatedData = zeros(numel(uniqReadID),(4*6) + 2);
    %change "6" to 5 in second argument if keeping the number of pairs
    %column in the new table
    oldTableData = table2array(currentTable(:,[2,3,6:9]));
    %iterate over the unique read IDs
    for i = 1:numel(uniqReadID)
        %match unique read IDs with the iteration number to generate a
        %reference index
        getUniqReadIDIdx = ic == i;
        %get repeated data
        getCurrentBins = getBinColumn(getUniqReadIDIdx);
        getCurrentWGBS = oldTableData(getUniqReadIDIdx, 1);
        getCurrentMethRatio = oldTableData(getUniqReadIDIdx, 2);
        %boolean checks to determine if the data was the same. python
        %script should always generate the same values for the WGBS and
        %meth ratio calculations within the same read ID.
        matchBinCheck = all(matches(getCurrentBins,binNames'));
        sameWGBSCheck = all(getCurrentWGBS(1) == getCurrentWGBS);
        sameMethRatioCheck = all(getCurrentMethRatio(1) == getCurrentMethRatio);
        %only include the values in the new table if the values were the
        %same for the given read ID
        if matchBinCheck && sameWGBSCheck && sameMethRatioCheck
            updatedData(i,:) = [getCurrentWGBS(1),...
                getCurrentMethRatio(1), ...
                reshape(oldTableData(getUniqReadIDIdx,3:end),1,[])];
        end
    end
    %fix variable names
    distDependentVarNames = reshape((getTableVariableNames(4:end)) + "_" + binNames',1,[]);
    %first three are "read_id","read_wgbs_distance", and
    %"read_meth_ratio". the remaining var names are the distance dependent
    %variables.
    updatedVarNames = [getTableVariableNames(1:3),distDependentVarNames];

    %remake the table
    updatedDataTable = array2table(updatedData,'VariableNames',updatedVarNames(2:end));
    rearrangedTable = [table(uniqReadID, 'VariableNames', updatedVarNames(1)), updatedDataTable];
end