% TRELIÇA PLANA NÃO LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% MATRIZ DE RIGIDEZ DO ELEMENTO PELA DEFORMAÇÃO DE ENGENHARIA
% 

function[Kc,Kg] = rigidez_engenharia_plana(Ae,Ee,Le,Lf,oxf,oyf,normal)

    % Matriz de rigidez constitutiva do elemento
    Kc = (Ee*Ae/Le)*[ oxf*oxf, oxf*oyf,-oxf*oxf,-oxf*oyf;
                      oyf*oxf, oyf*oyf,-oyf*oxf,-oyf*oyf;
                     -oxf*oxf,-oxf*oyf, oxf*oxf, oxf*oyf;
                     -oyf*oxf,-oyf*oyf, oyf*oxf, oyf*oyf];

    % Matriz de rigidez geométrica do elemento
    Kg = (normal/Lf)*[ oyf*oyf,     -oxf*oyf,     -oyf*oyf,       oxf*oyf;
                      -oyf*oxf,      oxf*oxf,      oyf*oxf,      -oxf*oxf;
                      -oyf*oyf,      oxf*oyf,      oyf*oyf,      -oxf*oyf;
                       oyf*oxf,     -oxf*oxf,     -oyf*oxf,       oxf*oxf];
                      
end