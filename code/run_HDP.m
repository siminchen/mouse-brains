
% Gamma priors for something...
alphaa = [1 1];
alphab = [1 1];
% Other parameters for Gibbs sampler
burn_in_samples = 1000;
num_samples = 100;
sample_spacing = 100;
num_initial_classes = 20;
% not sure what this one is yet
trainconparam = 1;

% Initializing hdpmix path stuff
initpath;


% Indicating multinomial distribution priors
func = hdpMultinomial_func;

% Reading data
disp('reading data')
d = dlmread('matlab_formatted_multinomialed_data.table');

% Useful to calculate once, and refer to it a lot
s = size(d);

disp('formatting/multinomializing...')
% Reformatting into cells
cell_d = cell(s(2),1);
for i = 1:length(cell_d)

	cell{i} = multinomialize_data(d(:,i)');
end

hyperparameter_priors = ones(s(1), 1) / s(1);
pause

disp('initializing HDP...')
hdp = hdp_init(func, 0,1, hyperparameter_priors, ...,
	alphaa, alphab);

[hdp dataindex] = hdp_adddp(hdp, length(cell_d) 1, 2);
hdp = hdp_setdata(hdp, dataindex, cell_d);

hdp = dp_activate(hdp, [1 dataindex], num_initial_classes)

disp('running sampling...')
[hdp sample lik] = hdp_posterior(hdp, burn_in_samples,...
	num_samples, sample_spacing, trainconparam, 1, 0, 1);

