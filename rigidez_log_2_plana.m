% TRELIÇA PLANA NÃO LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% MATRIZ DE RIGIDEZ DO ELEMENTO PELA DEFORMAÇÃO LOGARÍTMICA
% COM MUDANÇA NO VOLUME

function[Kc,Kg] = rigidez_log_2_plana(Ae,Ee,Le,Lf,oxf,oyf,normal,poisson)

    e_l = log(Lf/Le);               % Deformação logaritmica do elemento;

    % Matriz de rigidez constitutiva do elemento
    Kc = (Ee*Ae*Le*(e_l*(1-2*poisson)+1)/(Lf^2))*[ oxf*oxf, oxf*oyf,-oxf*oxf,-oxf*oyf;
                                                   oyf*oxf, oyf*oyf,-oyf*oxf,-oyf*oyf;
                                                  -oxf*oxf,-oxf*oyf, oxf*oxf, oxf*oyf;
                                                  -oyf*oxf,-oyf*oyf, oyf*oxf, oyf*oyf];

    % Matriz de rigidez geométrica do elemento
    Kg = (normal/Lf*(Le/Lf)^(2*poisson))*[ 1-2*oxf*oxf,   -2*oyf*oxf, -1+2*oxf*oxf,    2*oyf*oxf;
                                            -2*oxf*oyf,  1-2*oyf*oyf,    2*oxf*oyf, -1+2*oyf*oyf;
                                          -1+2*oxf*oxf,    2*oyf*oxf,  1-2*oxf*oxf,   -2*oyf*oxf;
                                             2*oxf*oyf, -1+2*oyf*oyf,   -2*oxf*oyf,  1-2*oyf*oyf];
                      
end