clear;

% set random seed 
rng(2)

% set paths
codeDir1 = '/gscratch/retina/vyomr/GitRepos/VLR-STRF/';  % YOUR CODE DIR HERE
codeDir2 = '/Users/riekelabbackup/Desktop/Vyom/gitrepos/VLR-STRF';  % YOUR CODE DIR HERE

dataPath = 'mea_data.mat';
% Load Y and X
mea_data = load(dataPath);
Y = mea_data.Y';
X = double(mea_data.X);
kSTA = double(mea_data.sta);
%%
try cd(codeDir1);
catch; cd(codeDir2);
end
set_paths;

% Data params
Nsamps = 19182; % number of time samples to in stimulus
%xdims = [95 152]; % spatial dimensions of stimulus. Only do red channel for now
xdims = [20 20];
nkx = prod(xdims);  % total number of spatial RF coeffs
nkt = 61;  % length of temporal filter (in bins)
dtbin = 1/120; % lenth of a single time bin
tmax = nkt*dtbin; % length of temporal RF

% build prior for spatial RF
spatPrior = build_vlrPrior('ALD',xdims);

% build prior for temporal RF
minlen_t = dtbin*2;   % minimum temporal lengthscale in normalised units
tempPrior = build_vlrPrior('TRD',nkt,minlen_t,tmax);
%kSTA = simpleRevcorr(X,Y,nkt);
% update initial hyperparameters from STA
[tempPrior, spatPrior] = initialiseHprs_vlrPriors(kSTA,tempPrior,spatPrior);

%% initialise model structure
rnk = 2;            % receptive field rank
opts = [];          % use default options

% build model structure 
m = build_vlrModel(Y,X,rnk,spatPrior,tempPrior,opts); 

%% Fit low-rank STRF using variational EM

% Set number of iterations per step of coordinate ascent
m.opts.maxiter.spatStep = 10;
m.opts.maxiter.tempStep = 10;
m.opts.maxiter.EM = 100;  % total number of EM iterations

% Run variational EM
fprintf('\nRunning variational EM...\n-------------------------\n\n');
m = fit_vlrModel(m);

%% Extract MAP filter estimate and hyper-parameter estimate

xhprs_hat = m.spatPrior.hprs;  % extract fitted spatial hyperparams
thprs_hat = m.tempPrior.hprs;  % extract fitted temporal hyperparams

% get maximum a posteriori estimate in first output argument
[mutHat,muxHat]  = getSTRF_vlrModel(m);

kMAP = mutHat*muxHat';
mut = mutHat*(mutHat\kt);  % representation of true components in temporal basis
mux = muxHat*(muxHat\kx);  % representation of true components in spatial basis

% Save the results
save('mea_202406_June.mat','kMAP','mut','mux','xhprs_hat','thprs_hat','m');