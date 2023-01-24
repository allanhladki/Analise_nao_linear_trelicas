% TRELIÇA PLANA NÃO LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% IMPLEMENTACAO ALGORITMO FINAL - TRELIÇAS PLANAS
%

clc
clear all %#ok<CLALL>

%% Dados iniciais

% Definir após os parâmetros necessários o nome do arquivo com as condições
% de entrada para análise;

[tipo_linear_fisico,tipo_deformacao,tipo_metodo_convergencia,nmax,imax,TOL,LAMBDA,DELTA_L,DELTA_L_0,DELTA_L_MAX,DELTA_L_MIN,comp_arco_variavel,Nd,zeta,plotar_deformada,plotar_trajetoria,plotar_t_d,GDL_trajetoria,Elemento_trajetoria,FE,plotcont,COORDNOS,CONEC,PROPELEM,REST,FC] = entrada();

%% MEF - Cálculos iniciais

NNOS = size(COORDNOS,1);                % Número de nós 
NBARRAS = size(PROPELEM,1);             % Número de barras
NNOSREST = size(REST,1);                % Número de nós restritos
GDL = 2*NNOS;                           % Número de graus de liberdade do problema
NELEM = size(CONEC,1);                  % Número de elementos
ESCOAMENTO = zeros(NELEM,1);            % Matriz do escoamento - Armazena quais barras estão escoando ou não
NNOSFC = size(FC,1);                    % Número de nós com forças concentradas
P = zeros(GDL,1);                       % Vetor de forças concentradas
K = zeros(GDL,GDL);                     % Matriz de rigidez global da estrutura
u = zeros(GDL,1);                       % Vetor dos deslocamentos globais da estrutura
normal = zeros(NELEM,1);                % Força atuando em cada elemento -> No início é zero
sigma = zeros(NELEM,1);                 % Tensão atuanda em cada elemento -> No início é zero
graph1 = zeros(nmax+1,5);               % Variável utilizada para geração de gráfico de resultados
SIGMAANT = zeros(NELEM);                % Tensão anterior do elemento - utilizado para verificar as deformações zero
fcc_i = zeros(GDL,1);                   % Vetor das cargas internas;
sigma_historia = zeros(NELEM,nmax+1);   % Vetor que armazena a historia das tensões de cada elemento;
e_historia = zeros(NELEM,nmax+1);       % Vetor que armazena a historia das deformações de cada elemento;
U=1;
        
% Alocar as FC no vetor P
for k = 1 : NNOSFC
    P(2*FC(k,1)-1) = FC(k,2);
    P(2*FC(k,1)) = FC(k,3);
end

if tipo_metodo_convergencia == 1 || tipo_metodo_convergencia == 2
    Delta_u = zeros(GDL,1);     % Vetor do incremento dos deslocamentos globais da estrutura

elseif tipo_metodo_convergencia == 3
    DELTA_U = zeros(GDL,1);             % Vetor do incremento dos deslocamentos globais da estrutura (SOMENTE PARA O MÉTODO DO COMP. DE ARCO);
    DELTA_U_1 = zeros(GDL,1);           % Vetor do incremento dos deslocamentos globais relacionado ao valor do fator de carga (COMP. ARCO HIPERPLANO FIXO)
    k_o = 0;                            % Alocação de variável para o cálculo do CST;
    DELTA_U_ANTERIOR = zeros(GDL,1);    % Vetor do acréscimo de deslocamentos do incremento anterior (SOMENTE PARA O MÉTODO DO COMP. DE ARCO);
    fcc_i_anterior = zeros(GDL,1);      % Vetor da forças internas do incremento anterior (SOMENTE PARA O MÉTODO DO COMP. DE ARCO);
    LAMBDA_ANTERIOR = 0.0;              % Alocação do LAMBDA anterior (SOMENTE PARA O MÉTODO DO COMP. DE ARCO);
end

% Plotar treliça indeformada
hold on
if plotar_deformada == 1
    for k = 1 : NELEM
        noi = CONEC(k,1);
        nof = CONEC(k,2);
        xi = COORDNOS(noi,1)+FE*u(2*noi-1);
        yi = COORDNOS(noi,2)+FE*u(2*noi);
        xf = COORDNOS(nof,1)+FE*u(2*nof-1);
        yf = COORDNOS(nof,2)+FE*u(2*nof);
        figure(1);
        line([xi,xf],[yi,yf],'Color','k','LineStyle','-','Marker','o','MarkerSize',5);
        title (sprintf('Configuração deformada da estrutura - Fator de escala = %i',FE));
    end
end

%% Processo Incremental-Iterativo

for n = 1:nmax

    if tipo_metodo_convergencia == 1
        Pp = (n/nmax)*P;
        for i = 1:imax
            [Kcc,fcc] = matriz_de_rigidez_plana(GDL,NELEM,CONEC,COORDNOS,u,PROPELEM,tipo_linear_fisico,tipo_deformacao,normal,Pp,NNOSREST,REST);
            [u,ACO_F,ACO_u,Delta_u] = newton_raphson_plano(fcc,Kcc,u,GDL,NELEM,CONEC,COORDNOS,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,Pp,NNOSREST,REST);
            [PROPELEM,sigma,normal,SIGMAANT,sigma_historia,e_historia] = esforcos_internos_plano(NELEM,CONEC,COORDNOS,u,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,n,sigma,normal,sigma_historia,e_historia);
            if ACO_F <= TOL && ACO_u <= TOL
                break;
            end
        end
        
        graph1(n+1,2) = (n/nmax);
        graph1(n+1,3) = u(GDL_trajetoria);
        
    elseif tipo_metodo_convergencia == 2
        Pp = (n/nmax)*P;
        [Kcc,fcc] = matriz_de_rigidez_plana(GDL,NELEM,CONEC,COORDNOS,u,PROPELEM,tipo_linear_fisico,tipo_deformacao,normal,Pp,NNOSREST,REST);
        for i = 1:imax
            [u,ACO_F,ACO_u,Delta_u] = newton_raphson_plano(fcc,Kcc,u,GDL,NELEM,CONEC,COORDNOS,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,Pp,NNOSREST,REST);
            
            if ACO_F <= TOL && ACO_u <= TOL
                break;
            end
        [PROPELEM,sigma,normal,SIGMAANT,sigma_historia,e_historia] = esforcos_internos_plano(NELEM,CONEC,COORDNOS,u,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,n,sigma,normal,sigma_historia,e_historia);
        end
        
        graph1(n+1,2) = (n/nmax);
        graph1(n+1,3) = u(GDL_trajetoria);
        
    elseif tipo_metodo_convergencia == 3
        Pp = P;
        [Kcc,fcc] = matriz_de_rigidez_plana(GDL,NELEM,CONEC,COORDNOS,u,PROPELEM,tipo_linear_fisico,tipo_deformacao,normal,Pp,NNOSREST,REST);
        [LAMBDA,u,DELTA_U,DELTA_L,k_o,CST,graph1,DELTA_U_ANTERIOR,fcc_i_anterior,LAMBDA_ANTERIOR] = comprimento_de_arco_NRM_plano(n,imax,TOL,Kcc,fcc,k_o,DELTA_L,DELTA_L_MAX,DELTA_L_MIN,DELTA_U,LAMBDA,u,comp_arco_variavel,Nd,zeta,GDL,NELEM,CONEC,COORDNOS,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,NNOSREST,REST,graph1,DELTA_U_ANTERIOR,fcc_i_anterior,LAMBDA_ANTERIOR);
        [PROPELEM,sigma,normal,SIGMAANT,sigma_historia,e_historia] = esforcos_internos_plano(NELEM,CONEC,COORDNOS,u,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,n,sigma,normal,sigma_historia,e_historia);
            
        graph1(n+1,1) = n;
        graph1(n+1,2) = LAMBDA;
        graph1(n+1,3) = u(GDL_trajetoria);
        graph1(n+1,5) = CST;
    end 
    
        % PLOTAGEM DA TRELIÇA DEFORMADA
    if mod(n,plotcont) == 0 || n == nmax
        if plotar_deformada == 1
            for k = 1 : NELEM
                noi = CONEC(k,1);
                nof = CONEC(k,2);
                xi = COORDNOS(noi,1)+FE*u(2*noi-1);
                yi = COORDNOS(noi,2)+FE*u(2*noi);
                xf = COORDNOS(nof,1)+FE*u(2*nof-1);
                yf = COORDNOS(nof,2)+FE*u(2*nof);
                figure(1);
                if n < nmax
                    line([xi,xf],[yi,yf],'Color','blue','LineStyle','--','Marker','o','MarkerSize',2);
                elseif n == nmax
                    line([xi,xf],[yi,yf],'Color','red','LineStyle','--','Marker','o','MarkerSize',2);
                end
            end
        end 
    end

end % Fim do processo incremental-iterativo


%% Saídas 

deslocamentos = [ [1:NNOS]' , vec2mat(u,2) ];
normal = [ [1:NELEM]' , normal ];
sigma = [ [1:NELEM]' , sigma ];

if plotar_trajetoria == 1
    figure(2);
    plot(-graph1(:,3),graph1(:,2));
    grid on;
    ylabel ('Fator de Carga','FontName','Cambria Math','FontSize',14);
    xlabel ('Deslocamento GDL 42','FontName','Cambria Math','FontSize',14);
end

if plotar_trajetoria == 1
    figure(3);
    plot(-graph1(:,3),graph1(:,5));
    grid on;
    ylabel ('CST','FontName','Cambria Math','FontSize',14);
    xlabel ('Deslocamento GDL 42','FontName','Cambria Math','FontSize',14);
end

if plotar_trajetoria == 1
    figure(4);
    plot(graph1(:,1),graph1(:,5));
    grid on;
    ylabel ('CST','FontName','Cambria Math','FontSize',14);
    xlabel ('Passos de Carga','FontName','Cambria Math','FontSize',14);
end


if plotar_t_d == 1
    figure(5); % Gráfico tensão x deformação
    plot(e_historia(Elemento_trajetoria,:),10*sigma_historia(Elemento_trajetoria,:));
    grid;
    title (sprintf('Curva tensão x deformação para o elemento %i',Elemento_trajetoria));
    xlabel ('Deformação');
    ylabel ('Tensão');
end
