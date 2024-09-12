global Lambda_w  Lambda_l  Lambda_r  Lambda_d  Lambda_tau  Lambda_tau_l ...
    C_w  C_l  C_r  C_d  C_tau  C_tau_l  A_H_w  A_H_l  A_H_r  A_H_d  A_H_tau  A_H_tau_l

if strcmp(model, '/rank') 

elseif strcmp(model, '/hank')

addpath([path vintage '/_auxiliary_functions/_hank_inputs'])

% get HANK Jacobians

if load_hank == 1
    load inputs_hank.mat
else
    get_hank_inputs_2
end

% make sure dimensions agree

Lambda_w     = Lambda_w(1:T,1:T);
Lambda_l     = Lambda_l(1:T,1:T);
Lambda_r     = Lambda_r(1:T,1:T);
Lambda_d     = Lambda_d(1:T,1:T);
Lambda_tau   = Lambda_tau(1:T,1:T);
Lambda_tau_l = Lambda_tau_l(1:T,1:T);

C_w     = C_w(1:T,1:T);
C_l     = C_l(1:T,1:T);
C_r     = C_r(1:T,1:T);
C_d     = C_d(1:T,1:T);
C_tau   = C_tau(1:T,1:T);
C_tau_l = C_tau_l(1:T,1:T);

A_H_w       = A_H_w(1:T,1:T);
A_H_l       = A_H_l(1:T,1:T);
A_H_r       = A_H_r(1:T,1:T);
A_H_d       = A_H_d(1:T,1:T);
A_H_tau     = A_H_tau(1:T,1:T);
A_H_tau_l   = A_H_tau_l(1:T,1:T);

else
    error('Select a valid model name')
end