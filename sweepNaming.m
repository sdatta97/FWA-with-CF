function naming = sweepNaming(sweepNames, sweepVals, aID)
%SWEEPNAMING automated result-naming convention.
%Every sweepable parameter is registered as a (name, values) pair:
%entries with numel > 1 are DETECTED as swept, scalars (or single-entry
%cells) as constant. The outputs encode:
%  folder : 'FWA_const_<names of the constant parameters>'
%  file   : 'results_<values of the constant parameters>_<seed>.csv'
%  header : CSV header = the SWEPT parameter names (their values fill
%           the table, one row per sweep-grid point) + record/outcome
%           columns
%  rowfmt : function handle mapping the registry-ordered CURRENT values
%           of all parameters to the swept-only CSV row prefix
%A truncated (verification) run turns axes into constants, so its
%outputs land in a DIFFERENT folder and can never mix with production.
assert(numel(sweepNames) == numel(sweepVals), 'registry length mismatch');
isSwept = cellfun(@(v) numel(v) > 1, sweepVals);
fmt1 = @(v) fmtVal(v);
naming.isSwept = isSwept;
naming.sweptNames = sweepNames(isSwept);
if any(~isSwept)
    naming.folder = ['FWA_const_' strjoin(sweepNames(~isSwept),'_')];
    cv = cellfun(fmt1, sweepVals(~isSwept), 'UniformOutput', false);
    naming.file = ['results_' strjoin(cv,'_') '_' aID '.csv'];
else
    naming.folder = 'FWA_const_none';
    naming.file = ['results_' aID '.csv'];
end
naming.header = [strjoin([sweepNames(isSwept), {'numCPE','numUE','init_FWA','Band_FWA','cell_se_ue','FWA_se'}], ',') '\n'];
naming.rowfmt = @(currVals) strjoin(cellfun(fmt1, currVals(isSwept), 'UniformOutput', false), ',');
end

function s = fmtVal(v)
%scalar parameter value -> filename/CSV token
if iscell(v), v = v{1}; end
if ischar(v) || isstring(v)
    s = char(v);
else
    s = num2str(v, '%g');
end
end
