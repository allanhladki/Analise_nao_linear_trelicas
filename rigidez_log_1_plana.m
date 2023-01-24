% TRELI�A PLANA N�O LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOL�GICA FEDERAL DO PARAN�
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% MATRIZ DE RIGIDEZ DO ELEMENTO PELA DEFORMA��O LOGAR�TMICA
% SEM MUDAN�A NO VOLUME

function[Kc,Kg] = rigidez_log_1_plana(Ae,Ee,Le,Lf,oxf,oyf,normal)

    % Matriz de rigidez constitutiva do elemento
    Kc = (Ee*Ae*Le/(Lf^2))*[ oxf*oxf, oxf*oyf,-oxf*oxf,-oxf*oyf;
                             oyf*oxf, oyf*oyf,-oyf*oxf,-oyf*oyf;
                            -oxf*oxf,-oxf*oyf, oxf*oxf, oxf*oyf;
                            -oyf*oxf,-oyf*oyf, oyf*oxf, oyf*oyf];

    % Matriz de rigidez geom�trica do elemento
    Kg = (normal*Le/Lf^2)*[ 1-2*oxf*oxf,   -2*oyf*oxf, -1+2*oxf*oxf,    2*oyf*oxf;
                             -2*oxf*oyf,  1-2*oyf*oyf,    2*oxf*oyf, -1+2*oyf*oyf;
                           -1+2*oxf*oxf,    2*oyf*oxf,  1-2*oxf*oxf,   -2*oyf*oxf;
                              2*oxf*oyf, -1+2*oyf*oyf,   -2*oxf*oyf,  1-2*oyf*oyf];
                      
end