%% FORECASTING WITH SW FACTORS
% Tomas Caravello, Alisdair McKay, Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

task = '/var_inputs';

%% SETTINGS

% VAR specification

n_lags     = 2;                    % number of lags
constant   = 0;                    % constant?
IRF_hor    = 20;

% dates for predictions

startdate = 1960.5;
enddate   = 2015;

%% FORECASTING RUNS

%----------------------------------------------------------------
% 1-Quarter-Ahead
%----------------------------------------------------------------

R2Table = zeros(3,2,4);

HORIZON = 1;  % this is the horizon of the forecast: y_t = A y_{t-HORIZON} + B y_{t-HORIZON-1} + C Factors_{t-HORIZON}
sw_factors_core

%----------------------------------------------------------------
% 4-Quarter-Ahead
%----------------------------------------------------------------

HORIZON = 4;  
sw_factors_core

%% FINAL TABLE

cd([path vintage task])

f = fopen(fullfile('_results', strcat('sw_compare', '.tex')), 'w'); % open file for writing

fprintf(f,'\\begin{table}\n');
fprintf(f,'\\centering\n');
fprintf(f,'\\begin{tabular}{lccccc} \\toprule \\toprule\n');
fprintf(f,'              & \\multicolumn{2}{c}{\\textbf{1-quarter ahead}}              &        & \\multicolumn{2}{c}{\\textbf{4-quarter ahead}}                         \\\\\n');
fprintf(f,'\\cline{2-3} \\cline{5-6}\n');
fprintf(f,'              & w/o F & w/ F & & w/o F & w/ F \\\\\n');
fprintf(f,'\\midrule\n');

% Define the data rows
data = {
    'Output gap', R2Table(1,1,1), R2Table(1,2,1), '', R2Table(1,1,4), R2Table(1,2,4);
    'Inflation', R2Table(2,1,1), R2Table(2,2,1), '', R2Table(2,1,4), R2Table(2,2,4);
    'Interest rate', R2Table(3,1,1), R2Table(3,2,1), '', R2Table(3,1,4), R2Table(3,2,4)
};

% Loop through each row of the data and print it
for i = 1:size(data, 1)
    fprintf(f,'%s & %.3f & %.3f & %s & %.3f & %.3f \\\\\n', data{i, :});
end

% Close the table
fprintf(f,'\\bottomrule\n');
fprintf(f,'\\end{tabular}\n');
fprintf(f,'\\caption{Assessing the incremental information content of the Stock-Watson factors: forecasting $R^2$ with and without inclusion of factors in the VAR.}\n');
fprintf(f,'\\label{tab:stock_watson_compare}\n');
fprintf(f,'\\end{table}\n');

fclose(f);