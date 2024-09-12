function sd = sd_fun(omega,Theta_Wold)

% preparations
hor = size(Theta_Wold,1);

% derive Theta(omega)
Theta_aux = 0;
for l = 1:hor
    Theta_aux = Theta_aux + squeeze(Theta_Wold(l,:,:)) * exp(-1i * omega * l);
end

% derive final result
sd = Theta_aux * ctranspose(Theta_aux);

end