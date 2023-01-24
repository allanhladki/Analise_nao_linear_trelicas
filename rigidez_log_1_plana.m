% TRELIÇA PLANA NÃO LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% MATRIZ DE RIGIDEZ DO ELEMENTO PELA DEFORMAÇÃO LOGARÍTMICA
% SEM MUDANÇA NO VOLUME

function[Kc,Kg] = rigidez_log_1_plana(Ae,Ee,Le,Lf,oxf,oyf,normal)

    % Matriz de rigidez constitutiva do elemento
    Kc = (Ee*Ae*Le/(Lf^2))*[ oxf*oxf, oxf*oyf,-oxf*oxf,-oxf*oyf;
                             oyf*oxf, oyf*oyf,-oyf*oxf,-oyf*oyf;
                            -oxf*oxf,-oxf*oyf, oxf*oxf, oxf*oyf;
                            -oyf*oxf,-oyf*oyf, oyf*oxf, oyf*oyf];

    % Matriz de rigidez geométrica do elemento
    Kg = (normal*Le/Lf^2)*[ 1-2*oxf*oxf,   -2*oyf*oxf, -1+2*oxf*oxf,    2*oyf*oxf;
                             -2*oxf*oyf,  1-2*oyf*oyf,    2*oxf*oyf, -1+2*oyf*oyf;
                           -1+2*oxf*oxf,    2*oyf*oxf,  1-2*oxf*oxf,   -2*oyf*oxf;
                              2*oxf*oyf, -1+2*oyf*oyf,   -2*oxf*oyf,  1-2*oyf*oyf];
                      
end