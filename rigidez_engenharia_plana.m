% TRELI�A PLANA N�O LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOL�GICA FEDERAL DO PARAN�
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% MATRIZ DE RIGIDEZ DO ELEMENTO PELA DEFORMA��O DE ENGENHARIA
% 

function[Kc,Kg] = rigidez_engenharia_plana(Ae,Ee,Le,Lf,oxf,oyf,normal)

    % Matriz de rigidez constitutiva do elemento
    Kc = (Ee*Ae/Le)*[ oxf*oxf, oxf*oyf,-oxf*oxf,-oxf*oyf;
                      oyf*oxf, oyf*oyf,-oyf*oxf,-oyf*oyf;
                     -oxf*oxf,-oxf*oyf, oxf*oxf, oxf*oyf;
                     -oyf*oxf,-oyf*oyf, oyf*oxf, oyf*oyf];

    % Matriz de rigidez geom�trica do elemento
    Kg = (normal/Lf)*[ oyf*oyf,     -oxf*oyf,     -oyf*oyf,       oxf*oyf;
                      -oyf*oxf,      oxf*oxf,      oyf*oxf,      -oxf*oxf;
                      -oyf*oyf,      oxf*oyf,      oyf*oyf,      -oxf*oyf;
                       oyf*oxf,     -oxf*oxf,     -oyf*oxf,       oxf*oxf];
                      
end