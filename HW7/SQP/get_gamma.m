function [gamma_k] = get_gamma(Lgrad1,Lgrad2)
    gamma_k = Lgrad2 - Lgrad1;
end