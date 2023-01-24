% TRELI�A PLANA N�O LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOL�GICA FEDERAL DO PARAN�
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% DADOS DE ENTRADA
% 

function[tipo_linear_fisico,tipo_deformacao,tipo_metodo_convergencia,nmax,imax,TOL,LAMBDA,DELTA_L,DELTA_L_0,DELTA_L_MAX,DELTA_L_MIN,comp_arco_variavel,Nd,zeta,plotar_deformada,plotar_trajetoria,plotar_t_d,GDL_trajetoria,Elemento_trajetoria,FE,plotcont,COORDNOS,CONEC,PROPELEM,REST,FC] = entrada()

% DEFINIR O TIPO DE AN�LISE F�SICA:
% 1 - An�lise linear sem limite de escoamento - Lei de Hooke;
% 2 - An�lise n�o linear f�sica - Modelo bilinear;
tipo_linear_fisico = 1;

% DEFINIR O TIPO DE DEFORMA��O UTILIZADA:
% 1 - Deforma��o de Engenharia;
% 2 - Deforma��o de Green;
% 3 - Deforma��o Logaritmica sem mudan�a de volume;
% 4 - Deforma��o Logaritmica com mudan�a de volume - OBS: Somente utilizar dentro do limite de escoamento;
tipo_deformacao = 1;

% DEFINIR O TIPO DE M�TODO DE CONVERG�NCIA
% 1 - M�todo de Newton-Raphson;
% 2 - M�todo de Newton-Raphson Modificado;
% 3 - M�todo do Comprimento de Arco de Riks-Wempner - NRM;
tipo_metodo_convergencia = 3;

% Dados de entrada para o m�todo incremetal-iterativo
nmax = 2750;               % N�mero m�ximo de incrementos; 
imax = 10;                 % N�mero m�ximo de itera��es;
TOL = 1e-6;                % Toler�ncia da converg�ncia;
LAMBDA = 0;                % Fator de carga (SOMENTO PARA O M�TODO DO COMP. DE ARCO);
DELTA_L = 0.1;             % Comprimento de arco inicial (SOMENTE PARA O M�TODO DO COMP. DE ARCO);     
DELTA_L_0 = DELTA_L;       % Comprimento de arco fixo para ajuste do comprimento de cada incremento (SOMENTE PARA O M�TODO DO COMP. DE ARCO);
DELTA_L_MAX = 1;           % Maior valor poss�vel para o comprimento de arco;
DELTA_L_MIN = 0.0001;      % Menor valor poss�vel para o comprimento de arco; 
comp_arco_variavel = 1;    % Ajuste do comprimento de arco (1 - DELTA_L VARI�VEL ����� 0 -> DELTA_L FIXO);
Nd = 3;                    % N�mero desejado de itera��es para cada ciclo (SOMENTE PARA M�TODO DE COMPRIMENTO DE ARCO);
zeta = 0.5;                % Fator da raz�o Nd/(Ni-1) (SOMENTE PARA M�TODO DE COMPRIMENTO DE ARCO);

% Dados de entrada para saida visual
plotar_deformada = 1;       % Plotagem da treli�a na configura��o deformada (1 -> PLOTAR ����� 0 -> N�O PLOTAR);
plotar_trajetoria = 1;      % Plotagem da trajet�ria de equilibrio para o GDL_trajetoria (1 -> PLOTAR ����� 0 -> N�O PLOTAR);
plotar_t_d = 1;             % Plotagem do grafico de tens�o x deforma��o para o Elemento_trajetoria (1 -> PLOTAR ����� 0 -> N�O PLOTAR);
GDL_trajetoria = 42;        % Grau de liberdade acompanhado na trajet�ria de equilibrio gr�fica;
Elemento_trajetoria = 1;    % Elemento acompanhado na trajetoria de tens�o x deforma��o;
FE = 1;                     % Fator de escala para plotagem da configura��o deformada;
plotcont = 1000;            % N�mero de incrementos entre cada plotagem da configura��o deformada;

% GEOMETRIA
% Coordenadas dos n�s - [coordX (cm), coordY (cm)]
COORDNOS = [-33.9411,	33.9411;
            -35.3553,	35.3553;
            -31.1735,	36.4995;
            -32.4724,	38.0203;
            -28.2137,	38.8328;
            -29.3893,	40.4508;
            -25.0799,	40.9267;
            -26.1249,	42.6320;
            -21.7915,	42.7683;
            -22.6995,	44.5503;
            -18.3688,	44.3462;
            -19.1342,	46.1940;
            -14.8328,	45.6507;
            -15.4508,	47.5528;
            -11.2054,	46.6738;
            -11.6723,	48.6185;
             -7.5089,	47.4090;
             -7.8217,	49.3844;
             -3.7660,	47.8520;
             -3.9230,	49.8459;
              0.0000,	48.0000;
              0.0000,	50.0000;
              3.7660,	47.8520;
              3.9230,	49.8459;
              7.5089,	47.4090;
              7.8217,	49.3844;
             11.2054,	46.6738;
             11.6723,	48.6185;
             14.8328,	45.6507;
             15.4508,	47.5528;
             18.3688,	44.3462;
             19.1342,	46.1940;
             21.7915,	42.7683;
             22.6995,	44.5503;
             25.0799,	40.9267;
             26.1249,	42.6320;
             28.2137,	38.8328;
             29.3893,	40.4508;
             31.1735,	36.4995;
             32.4724,	38.0203;
             33.9411,	33.9411;
             35.3553,	35.3553;];

% Barras - Matriz de Conectividade [no inicial, n� final~]
CONEC = [   1,  3;
            3,  5;
            5,  7;
            7,  9;
            9, 11;
           11, 13;
           13, 15;
           15, 17;
           17, 19;
           19, 21;
           21, 23;
           23, 25;
           25, 27;
           27, 29;
           29, 31;
           31, 33;
           33, 35;
           35, 37;
           37, 39;
           39, 41;
            2,  4;
            4,  6;
            6,  8;
            8, 10;
           10, 12;
           12, 14;
           14, 16;
           16, 18;
           18, 20;
           20, 22;
           22, 24;
           24, 26;
           26, 28;
           28, 30;
           30, 32;
           32, 34;
           34, 36;
           36, 38;
           38, 40;
           40, 42;
            1,  4;
            3,  6;
            5,  8;
            7, 10;
            9, 12;
           11, 14;
           13, 16;
           15, 18;
           17, 20;
           19, 22;
           21, 24;
           23, 26;
           25, 28;
           27, 30;
           29, 32;
           31, 34;
           33, 36;
           35, 38;
           37, 40;
           39, 42;
            2,  3;
            4,  5;
            6,  7;
            8,  9;
           10, 11;
           12, 13;
           14, 15;
           16, 17;
           18, 19;
           20, 21;
           22, 23;
           24, 25;
           26, 27;
           28, 29;
           30, 31;
           32, 33;
           34, 35;
           36, 37;
           38, 39;
           40, 41;
            1,  2;
            3,  4;
            5,  6;
            7,  8;
            9, 10;
           11, 12;
           13, 14;
           15, 16;
           17, 18;
           19, 20;
           21, 22;
           23, 24;
           25, 26;
           27, 28;
           29, 30;
           31, 32;
           33, 34;
           35, 36;
           37, 38;
           39, 40;
           41, 42;];
           
% Se��o transversal
A1 = 100;          % cm�
 
% Material
E1 = 5e5;          % kN/cm^2
v = 0.3;

% M�dulo de elasticidade ap�s o escoamento do material
% Somente necess�rio para an�lise n�o linear f�sica
H = 0;           % kN/cm^2
Et = (H*E1)/(H+E1);

% Tens�o de escoamento do material para an�lise n�o linear f�sica
% A partir dessa tens�o o material estar� em regime pl�stico
% Somente necess�rio para an�lise n�o linear f�sica
TES1 = 15;          % kN/cm�

% Barras - Propriedades
% �rea da se��o;
% M�dulo de elasticidade el�stico;
% M�dulo de elasticidade pl�stico (SOMENTE PARA AN�LISE N�O LINEAR F�SICA);
% Tens�o de escoamento (SOMENTE PARA AN�LISE N�O LINEAR F�SICA);
% Deforma��o associada a tens�o de escoamento (SOMENTE PARA AN�LISE N�O LINEAR F�SICA);
% e_p, deforma��o zero - a �ltima deforma��o tal que a tens�o seja 0 (SOMENTE PARA AN�LISE N�O LINEAR F�SICA);
% Contador para c�lculo das tens�es de escoamento -> Sempre igual a 0 (SOMENTE PARA AN�LISE N�O LINEAR F�SICA);
% Coeficiente de Poisson (SOMENTE PARA MODELO COM MUDAN�A DE VOLUME);
PROPELEM = [A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;
            A1,E1,Et,TES1,TES1/E1,0,0,v;];
                
% Restri��es - [N�, RESTX(0/1), RESTY(0/1)]
REST = [1, 1,1;
        41,1,1;];
   
% For�as externas - [N�, FX(kN), FY(kN)]
FC = [22,0,-1];

end
