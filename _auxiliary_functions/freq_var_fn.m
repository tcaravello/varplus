function freq_var = freq_var_fn(Theta,omega_1,omega_2);

IRF_hor = size(Theta,3);

n_omega = 1e2;
omega_grid = linspace(omega_1,omega_2,n_omega)';
domega = omega_grid(2) - omega_grid(1);
freq_var = 0;
for i_omega = 1:n_omega
    temp = 0;
    for i_hor = 1:IRF_hor
        temp = temp + Theta(:,:,i_hor) * exp(-1i*omega_grid(i_omega)*(i_hor-1));
    end
    freq_var = freq_var + 1/(2*pi) * temp * ctranspose(temp) * domega;
end
freq_var = real(freq_var);