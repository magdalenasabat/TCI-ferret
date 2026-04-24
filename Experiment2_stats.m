% Loads published data and reproduces correlation computation and model
% fitting.
% 
% Sabat, M., Gouyette, H., Gaucher, Q., Espejo, M. L., David, S. V., Norman-Haignere, S., & Boubenec, Y. (2025). 
% Neurons in auditory cortex integrate information within a constrained and context-invariant temporal window. 
% Current Biology, 0(0). https://doi.org/10.1016/j.cub.2025.11.011
% 
%
% author: Magdalena email: magdalena.sabat@ens.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize
clc; close all; clear all; 


addpath(genpath(strcat(['code' filesep 'myFunctions'])));
addpath(genpath(strcat(['code' filesep 'toolbox-sam-tci'])));

dataFolder = strcat([ 'data' filesep 'experiment2']);

subnames={'D','E','F','G','H'};
sessions = {[17:29,31,33], ...
    [29:32],...
    [3:5,11,12,18,19],...
    [9,10,15:18,23,24],...
    [10,21,23]};

% for plotting 
markerstyles = {'^','square', 'o'}; 

%% Get all data

xi=1;
all_widths = zeros(1,500);
all_centers = zeros(1,500);
all_layers = zeros(1,500);
all_areas = cell(1,500);
all_subs = zeros(1,500);
all_unique_sess = cell(1,500);
 


for subid= 1:length(subnames) 
    
    subjid = subnames{subid};
    fprintf(strcat([' SUBJECT ' subjid '\n']))
    
    
    for sessN=sessions{subid}
        
%         fprintf(strcat([' sess #' sprintf('%03d',sessN) '\n']))

        %load files
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_anatomical.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_model_fits.mat'  ]));
        load(strcat([ dataFolder filesep 'sub-' subjid filesep 'sess-' sprintf('%03d',sessN) '_cross_context_correlation.mat'  ]));
 

        % pickup model fits and distance
        units_n = length(M.best_intper_sec);
        all_widths(xi:xi+units_n-1)=M.best_intper_sec;
        all_centers(xi:xi+units_n-1)=M.best_delay_sec_median;
        all_subs(xi:xi+units_n-1)=subid;
        all_layers(xi:xi+units_n-1)=layer;
        all_areas(xi:xi+units_n-1)=area;
        
        
        all_unique_sess(xi:xi+units_n-1) = repmat( {[sprintf('%03d',sessN) subjid ]}, 1, units_n);
        xi=xi+units_n;
        
    end
end


all_widths = all_widths(1:xi-1); all_widths=all_widths*1000;
all_centers = all_centers(1:xi-1); all_centers=all_centers*1000;
all_layers = all_layers(1:xi-1);
all_areas = all_areas(1:xi-1);
all_subs = all_subs(1:xi-1);
all_unique_sess = all_unique_sess(1:xi-1);

%%
% turn layers to string
layers_string = arrayfun(@(x) num2str(x), all_layers , 'UniformOutput', false);

%create the table
tbl = table(log2(all_widths(:)),log2(all_centers(:)),all_areas(:),layers_string(:),all_unique_sess(:),'VariableNames',{'width','center','region','layer','recording'});

%run the lme
lme = fitlme(tbl,'width ~ 1 + region + layer + (1| recording)',  'CovariancePattern', 'diagonal')
% test signigicance of fixed effects
anova(lme,'DFMethod', 'satterthwaite')

%%
% binarize layer variable
layers_binary = zeros(size(all_layers));
layers_binary(all_layers==4)=1;
layers_binary(all_layers==13)=1;
% turn into categorical, the second argument changes the direction of the
% comparison in the lme
layers_binary=categorical(layers_binary,[0 1]);

%create the table
tbl = table(log2(all_widths(:)),log2(all_centers(:)),all_areas(:),layers_binary(:),all_unique_sess(:),'VariableNames',{'width','center','region','layer','recording'});

%run the lme
lme = fitlme(tbl,'width ~ 1 + region + layer + (1| recording)',  'CovariancePattern', 'diagonal');



% test the significance of the difference of the fixed effects of layer and
% region

H = [0, 1, -1]; % picks out region and layer 
C = [0];
[p1, F, df1, df2] = coefTest(lme, H, C, 'DFMethod', 'satterthwaite');
disp(['Significance of the difference between the region and layer effects: FStat = ' num2str(F) ', DF1 = ' num2str(df1) ', DF2 = ' num2str(df2)  ',p = ' num2str(p1)] )
H = [0, 1, 1]; % picks out region and layer 
C = [0];
[p2, F, df1, df2] = coefTest(lme, H, C, 'DFMethod', 'satterthwaite');
disp(['Significance of the sum between the region and layer effects: FStat = ' num2str(F) ', DF1 = ' num2str(df1) ', DF2 = ' num2str(df2)  ',p = ' num2str(p2)] )
disp(['The union of the two hypotheses, max(p) = ' num2str(max(p1,p2))])
disp(['Absolute difference between the layers when grouped, M = ' num2str( median(all_widths(layers_binary==categorical(1)))-median(all_widths(layers_binary==categorical(-1))) ) 'ms'] )


