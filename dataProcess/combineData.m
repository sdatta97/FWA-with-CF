%Taken from 
% https://in.mathworks.com/matlabcentral/answers/
%   538119-how-to-import-to-matlab-many-text-files-as-table-type
projectdir = '/Users/sdatta/FWA-with-CF/resultData/ratecomparisonResults';
dinfo = dir(fullfile(projectdir, '*.csv'));   %use appropriate extension
filenames = fullfile({dinfo.folder}, {dinfo.name});
nfiles = length(filenames);
tables = cell(nfiles,1);
for K = 1 : nfiles
    tables{K} = readtable(filenames{K});
end

combinedTable = vertcat(tables{:});

colNames = combinedTable.Properties.VariableNames;
%% group the same parameter iterations and get stats mean and std, can add other stuff
% changingVars = cell(1,length(colNames)-6);
% for i=1:(length(colNames)-6)
%     changingVars{i} = colNames{i};
% end
% changingVars = cell(1,length(colNames)-4);
% for i=1:(length(colNames)-4)
%     changingVars{i} = colNames{i};
% end
changingVars = cell(1,length(colNames)-5);
for i=1:(length(colNames)-5)
    changingVars{i} = colNames{i};
end
% changingVars = cell(1,length(colNames)-2);
% for i=1:(length(colNames)-2)
%     changingVars{i} = colNames{i};
% end

summaryTable  = groupsummary(combinedTable,changingVars,{'mean','std','median'});

writetable(summaryTable,'./rate_comp_data_vary_pf_bs_density_5.txt')
writetable(summaryTable,'./rate_comp_data_vary_pf_bs_density_5.csv')