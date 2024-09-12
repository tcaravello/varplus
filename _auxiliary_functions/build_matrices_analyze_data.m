function [Pi_m_struct,Y_m_struct,R_n_m_struct, m_fit_struct,Pi_m_fit_struct,Y_m_fit_struct, R_n_m_fit_struct] =...
    build_matrices_analyze_data(Pi_m_collector, Y_m_collector, R_n_m_collector, m_fit_collector, Pi_m_fit_collector,Y_m_fit_collector, R_n_m_fit_collector)

T_save = size(Pi_m_collector,1);
T_use   = size(m_fit_collector,1);
Pi_m_struct = struct();
Y_m_struct = struct();
R_n_m_struct = struct();
m_fit_struct = struct();

% percentiles 0.05; 16; 0.5; 0.84; 0.95
Pi_m_struct.lb = zeros(T_save,T_save);
Pi_m_struct.ub = zeros(T_save,T_save);
Pi_m_struct.median = zeros(T_save,T_save);
Pi_m_struct.per_05 = zeros(T_save,T_save);
Pi_m_struct.per_95 = zeros(T_save,T_save);
Pi_m_struct.mean = zeros(T_save,T_save);
Pi_m_struct.sd = zeros(T_save,T_save);


Y_m_struct.lb = zeros(T_save,T_save);
Y_m_struct.ub = zeros(T_save,T_save);
Y_m_struct.median = zeros(T_save,T_save);
Y_m_struct.per_05 = zeros(T_save,T_save);
Y_m_struct.per_95 = zeros(T_save,T_save);
Y_m_struct.mean = zeros(T_save,T_save);
Y_m_struct.sd = zeros(T_save,T_save);

R_n_m_struct.lb = zeros(T_save,T_save);
R_n_m_struct.ub = zeros(T_save,T_save);
R_n_m_struct.median = zeros(T_save,T_save);
R_n_m_struct.per_05 = zeros(T_save,T_save);
R_n_m_struct.per_95 = zeros(T_save,T_save);
R_n_m_struct.mean = zeros(T_save,T_save);
R_n_m_struct.sd = zeros(T_save,T_save);

m_fit_struct.lb = zeros(T_use,1);
m_fit_struct.ub = zeros(T_use,1);
m_fit_struct.median = zeros(T_use,1);
m_fit_struct.per_05 = zeros(T_use,1);
m_fit_struct.per_95 = zeros(T_use,1);
m_fit_struct.mean = zeros(T_use,1);
m_fit_struct.sd = zeros(T_use,1);

% fitted IRFs
Pi_m_fit_struct.lb = zeros(T_use,1);
Pi_m_fit_struct.ub = zeros(T_use,1);
Pi_m_fit_struct.median = zeros(T_use,1);
Pi_m_fit_struct.per_05 = zeros(T_use,1);
Pi_m_fit_struct.per_95 = zeros(T_use,1);
Pi_m_fit_struct.mean = zeros(T_use,1);
Pi_m_fit_struct.sd = zeros(T_use,1);


Y_m_fit_struct.lb = zeros(T_use,1);
Y_m_fit_struct.ub = zeros(T_use,1);
Y_m_fit_struct.median = zeros(T_use,1);
Y_m_fit_struct.per_05 = zeros(T_use,1);
Y_m_fit_struct.per_95 = zeros(T_use,1);
Y_m_fit_struct.mean = zeros(T_use,1);
Y_m_fit_struct.sd = zeros(T_use,1);

R_n_m_fit_struct.lb = zeros(T_use,1);
R_n_m_fit_struct.ub = zeros(T_use,1);
R_n_m_fit_struct.median = zeros(T_use,1);
R_n_m_fit_struct.per_05 = zeros(T_use,1);
R_n_m_fit_struct.per_95 = zeros(T_use,1);
R_n_m_fit_struct.mean = zeros(T_use,1);
R_n_m_fit_struct.sd = zeros(T_use,1);




for ii=1:T_use
            m_fit_struct.median(ii) = quantile(m_fit_collector(ii,:), 0.5);
            m_fit_struct.ub(ii) = quantile(m_fit_collector(ii,:),0.84);
            m_fit_struct.lb(ii) = quantile(m_fit_collector(ii,:),0.16);
            m_fit_struct.per_05(ii) = quantile(m_fit_collector(ii,:),0.05);
            m_fit_struct.per_95(ii) = quantile(m_fit_collector(ii,:),0.95);
            m_fit_struct.mean(ii) = mean(m_fit_collector(ii,:));
            m_fit_struct.sd(ii) = std(m_fit_collector(ii,:));

            Pi_m_fit_struct.median(ii) = quantile(Pi_m_fit_collector(ii,:), 0.5);
            Pi_m_fit_struct.ub(ii) = quantile(Pi_m_fit_collector(ii,:),0.84);
            Pi_m_fit_struct.lb(ii) = quantile(Pi_m_fit_collector(ii,:),0.16);
            Pi_m_fit_struct.per_05(ii) = quantile(Pi_m_fit_collector(ii,:),0.05);
            Pi_m_fit_struct.per_95(ii) = quantile(Pi_m_fit_collector(ii,:),0.95);
            Pi_m_fit_struct.mean(ii) = mean(Pi_m_fit_collector(ii,:));
            Pi_m_fit_struct.sd(ii) = std(Pi_m_fit_collector(ii,:));

            Y_m_fit_struct.median(ii) = quantile(Y_m_fit_collector(ii,:), 0.5);
            Y_m_fit_struct.ub(ii) = quantile(Y_m_fit_collector(ii,:),0.84);
            Y_m_fit_struct.lb(ii) = quantile(Y_m_fit_collector(ii,:),0.16);
            Y_m_fit_struct.per_05(ii) = quantile(Y_m_fit_collector(ii,:),0.05);
            Y_m_fit_struct.per_95(ii) = quantile(Y_m_fit_collector(ii,:),0.95);
            Y_m_fit_struct.mean(ii) = mean(Y_m_fit_collector(ii,:));
            Y_m_fit_struct.sd(ii) = std(Y_m_fit_collector(ii,:));

            R_n_m_fit_struct.median(ii) = quantile(R_n_m_fit_collector(ii,:), 0.5);
            R_n_m_fit_struct.ub(ii) = quantile(R_n_m_fit_collector(ii,:),0.84);
            R_n_m_fit_struct.lb(ii) = quantile(R_n_m_fit_collector(ii,:),0.16);
            R_n_m_fit_struct.per_05(ii) = quantile(R_n_m_fit_collector(ii,:),0.05);
            R_n_m_fit_struct.per_95(ii) = quantile(R_n_m_fit_collector(ii,:),0.95);
            R_n_m_fit_struct.mean(ii) = mean(R_n_m_fit_collector(ii,:));
            R_n_m_fit_struct.sd(ii) = std(R_n_m_fit_collector(ii,:));

    for jj=1:T_save
            Pi_m_struct.median(ii,jj) = quantile(Pi_m_collector(ii,jj,:), 0.5);
            Pi_m_struct.ub(ii,jj) = quantile(Pi_m_collector(ii,jj,:),0.84);
            Pi_m_struct.lb(ii,jj) = quantile(Pi_m_collector(ii,jj,:),0.16);
            Pi_m_struct.per_05(ii,jj) = quantile(Pi_m_collector(ii,jj,:),0.05);
            Pi_m_struct.per_95(ii,jj) = quantile(Pi_m_collector(ii,jj,:),0.95);
            Pi_m_struct.mean(ii,jj) = mean(Pi_m_collector(ii,jj,:));
            Pi_m_struct.sd(ii,jj) = std(Pi_m_collector(ii,jj,:));

            Y_m_struct.median(ii,jj) = quantile(Y_m_collector(ii,jj,:), 0.5);
            Y_m_struct.ub(ii,jj) = quantile(Y_m_collector(ii,jj,:),0.84);
            Y_m_struct.lb(ii,jj) = quantile(Y_m_collector(ii,jj,:),0.16);
            Y_m_struct.per_05(ii,jj) = quantile(Y_m_collector(ii,jj,:),0.05);
            Y_m_struct.per_95(ii,jj) = quantile(Y_m_collector(ii,jj,:),0.95);
            Y_m_struct.mean(ii,jj) = mean(Y_m_collector(ii,jj,:));
            Y_m_struct.sd(ii,jj) = std(Y_m_collector(ii,jj,:));


            R_n_m_struct.median(ii,jj) = quantile(R_n_m_collector(ii,jj,:), 0.5);
            R_n_m_struct.ub(ii,jj) = quantile(R_n_m_collector(ii,jj,:),0.84);
            R_n_m_struct.lb(ii,jj) = quantile(R_n_m_collector(ii,jj,:),0.16);
            R_n_m_struct.per_05(ii,jj) = quantile(R_n_m_collector(ii,jj,:),0.05);
            R_n_m_struct.per_95(ii,jj) = quantile(R_n_m_collector(ii,jj,:),0.95);
            R_n_m_struct.mean(ii,jj) = mean(R_n_m_collector(ii,jj,:));
            R_n_m_struct.sd(ii,jj) = std(R_n_m_collector(ii,jj,:));

    end
end
end