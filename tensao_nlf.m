% TRELIÇA PLANA NÃO LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% TENSÃO - MODELO BILINEAR 

function [sigma_atual,sigma_esc,e_esc,e_p,sigma_anterior,contador] = tensao_nlf(sigma_esc,e_esc,e_p,sigma_anterior,contador,E,Et,e_atual,i)
    
    sigma_prev = E*(e_atual-e_p);           % Verifica qual seria a tensão em regime linear
    
    if i==1 && sigma_prev < 0               % Muda o sinal da deformação de escoamento caso seja o primeiro passo de carga
        e_esc = -e_esc;
    elseif i ~= 1 && ((sigma_prev <0 && sigma_anterior >= 0) || (sigma_prev >=0 && sigma_anterior <0))
        e_esc = 2*e_p-e_esc; % Acha a deformação de escoamento com sinal contrário caso haja mudança no sentido da tensão
    end

    if sigma_prev > 0 && e_atual >= e_esc       % Ocorre escoamento para tensão positiva
            sigma_atual = sigma_esc + Et*(e_atual-e_esc);
            contador = 1;
            sigma_esc = sigma_atual;
            e_esc = e_atual;
    elseif sigma_prev < 0 && e_atual <= e_esc   % Ocorre escoamento para tensão negativa
            sigma_atual = -sigma_esc + Et*(e_atual-e_esc);
            contador = -1;
            sigma_esc = -sigma_atual;
            e_esc = e_atual;
    else         % Não ocorre escoamento
        if sigma_prev > 0
            sigma_atual = sigma_esc - E*(e_esc-e_atual);
        elseif sigma_prev == 0
            sigma_atual = 0;
        else
            sigma_atual = -sigma_esc - E*(e_esc-e_atual);
        end
        if contador == 1   % Atualiza o valor de e_p
            e_p = (sigma_esc*e_atual-e_esc*sigma_atual)/(sigma_esc-sigma_atual);
            contador = 0;
        elseif contador == -1
            e_p = (-sigma_esc*e_atual-e_esc*sigma_atual)/(-sigma_esc-sigma_atual);
            contador = 0;
        end
    end    
    
    sigma_anterior = sigma_atual;
    
end
