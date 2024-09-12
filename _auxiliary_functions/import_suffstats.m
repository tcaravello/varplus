%----------------------------------------------------------------
% RE Models
%----------------------------------------------------------------

if indic_RE == 1

% posterior across models

try
   load non_behav_models_draws
catch ME
   if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
      msg = ['The posterior draws for the policy shock matrices were not found. Please check that you correctly downloaded them from the Dropbox link in the Readme.'];
        causeException = MException('MATLAB:myCode:folder_not_found',msg);
        ME = addCause(ME,causeException);
   end
   rethrow(ME)
end



Pi_m_all = mean(4 * Pi_m_collector_non_behav,3);
I_m_all  = mean(4 * R_n_m_collector_non_behav,3);
Y_m_all  = mean(Y_m_collector_non_behav,3);

Pi_m_draws = 4 * Pi_m_collector_non_behav;
I_m_draws  = 4 * R_n_m_collector_non_behav;
Y_m_draws  = Y_m_collector_non_behav;

Pi_m_base  = mean(Pi_m_draws,3);
I_m_base   = mean(I_m_draws,3);
Y_m_base   = mean(Y_m_draws,3);

clear m_fit_collector_non_behav model_posterior_non_behav Pi_m_collector_non_behav R_n_m_collector_non_behav Y_m_collector_non_behav

% HANK

load hank_draws_main

Pi_m_hank = mean(4 * Pi_m_collector,3);
I_m_hank  = mean(4 * R_n_m_collector,3);
Y_m_hank  = mean(Y_m_collector,3);

Pi_m_hank_draws = 4 * Pi_m_collector;
I_m_hank_draws  = 4 * R_n_m_collector;
Y_m_hank_draws  = Y_m_collector;

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector

% RANK

load rank_draws_main

Pi_m_rank = mean(4 * Pi_m_collector,3);
I_m_rank  = mean(4 * R_n_m_collector,3);
Y_m_rank  = mean(Y_m_collector,3);

Pi_m_rank_draws = 4 * Pi_m_collector;
I_m_rank_draws  = 4 * R_n_m_collector;
Y_m_rank_draws  = Y_m_collector;

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector

end

%----------------------------------------------------------------
% Fixed Behavioral Models
%----------------------------------------------------------------

if indic_behav == 1

% posterior across models

try
   load behav_all_models_draws
catch ME
   if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
      msg = ['The posterior draws for the policy shock matrices were not found. Please check that you correctly downloaded them from the Dropbox link in the Readme.'];
        causeException = MException('MATLAB:myCode:folder_not_found',msg);
        ME = addCause(ME,causeException);
   end
   rethrow(ME)
end



Pi_m_all = mean(4 * Pi_m_collector_behav,3);
I_m_all  = mean(4 * R_n_m_collector_behav,3);
Y_m_all  = mean(Y_m_collector_behav,3);

Pi_m_draws = 4 * Pi_m_collector_behav;
I_m_draws  = 4 * R_n_m_collector_behav;
Y_m_draws  = Y_m_collector_behav;

Pi_m_base  = mean(Pi_m_draws,3);
I_m_base   = mean(I_m_draws,3);
Y_m_base   = mean(Y_m_draws,3);

clear m_fit_collector_behav model_posterior_behav Pi_m_collector_behav R_n_m_collector_behav Y_m_collector_behav

% HANK

load hank_draws_main_behav

Pi_m_hank = mean(4 * Pi_m_collector,3);
I_m_hank  = mean(4 * R_n_m_collector,3);
Y_m_hank  = mean(Y_m_collector,3);

Pi_m_hank_draws = 4 * Pi_m_collector;
I_m_hank_draws  = 4 * R_n_m_collector;
Y_m_hank_draws  = Y_m_collector;

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector

% RANK

load rank_draws_main_behav

Pi_m_rank = mean(4 * Pi_m_collector,3);
I_m_rank  = mean(4 * R_n_m_collector,3);
Y_m_rank  = mean(Y_m_collector,3);

Pi_m_rank_draws = 4 * Pi_m_collector;
I_m_rank_draws  = 4 * R_n_m_collector;
Y_m_rank_draws  = Y_m_collector;

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector
end

%----------------------------------------------------------------
% All Models Jointly
%----------------------------------------------------------------

if indic_joint == 1

% posterior across models

try
   load all_models_draws
catch ME
   if (strcmp(ME.identifier,'MATLAB:load:couldNotReadFile'))
      msg = ['The posterior draws for the policy shock matrices were not found. Please check that you correctly downloaded them from the Dropbox link in the Readme.'];
        causeException = MException('MATLAB:myCode:folder_not_found',msg);
        ME = addCause(ME,causeException);
   end
   rethrow(ME)
end




Pi_m_all = mean(4 * Pi_m_collector,3);
I_m_all  = mean(4 * R_n_m_collector,3);
Y_m_all  = mean(Y_m_collector,3);

Pi_m_draws = 4 * Pi_m_collector;
I_m_draws  = 4 * R_n_m_collector;
Y_m_draws  = Y_m_collector;

Pi_m_base  = mean(Pi_m_draws,3);
I_m_base   = mean(I_m_draws,3);
Y_m_base   = mean(Y_m_draws,3);

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector

% HANK

load hank_behav_and_non_behav_draws

Pi_m_hank = mean(4 * Pi_m_collector_only_hank,3);
I_m_hank  = mean(4 * R_n_m_collector_only_hank,3);
Y_m_hank  = mean(Y_m_collector_only_hank,3);

Pi_m_hank_draws = 4 * Pi_m_collector_only_hank;
I_m_hank_draws  = 4 * R_n_m_collector_only_hank;
Y_m_hank_draws  = Y_m_collector_only_hank;

clear m_fit_collector_only_hank model_posterior_only_hank Pi_m_collector_only_hank R_n_m_collector_only_hank Y_m_collector_only_hank

% RANK

load rank_behav_and_non_behav_draws

Pi_m_rank = mean(4 * Pi_m_collector_only_rank,3);
I_m_rank  = mean(4 * R_n_m_collector_only_rank,3);
Y_m_rank  = mean(Y_m_collector_only_rank,3);

Pi_m_rank_draws = 4 * Pi_m_collector_only_rank;
I_m_rank_draws  = 4 * R_n_m_collector_only_rank;
Y_m_rank_draws  = Y_m_collector_only_rank;

clear m_fit_collector_only_rank model_posterior_only_rank Pi_m_collector_only_rank R_n_m_collector_only_rank Y_m_collector_only_rank

% B-HANK

load hank_draws_main_behav

Pi_m_bhank = mean(4 * Pi_m_collector,3);
I_m_bhank  = mean(4 * R_n_m_collector,3);
Y_m_bhank  = mean(Y_m_collector,3);

Pi_m_bhank_draws = 4 * Pi_m_collector;
I_m_bhank_draws  = 4 * R_n_m_collector;
Y_m_bhank_draws  = Y_m_collector;

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector

% B-RANK

load rank_draws_main_behav

Pi_m_brank = mean(4 * Pi_m_collector,3);
I_m_brank  = mean(4 * R_n_m_collector,3);
Y_m_brank  = mean(Y_m_collector,3);

Pi_m_brank_draws = 4 * Pi_m_collector;
I_m_brank_draws  = 4 * R_n_m_collector;
Y_m_brank_draws  = Y_m_collector;

clear m_fit_collector model_posterior Pi_m_collector R_n_m_collector Y_m_collector

end