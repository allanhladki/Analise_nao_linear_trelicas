% TRELIÇA PLANA NÃO LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% MÉTODO DO COMPRIMENTO DE ARCO DO RIKS-WEMPNER
% CONTEMPLANDO ESTRATÉGIA DE CORREÇÃO

function[LAMBDA,u,DELTA_U,DELTA_L,k_o,CST,graph1,DELTA_U_ANTERIOR,fcc_i_anterior,LAMBDA_ANTERIOR] = comprimento_de_arco_NRM_plano(n,imax,TOL,Kcc,fcc,k_o,DELTA_L,DELTA_L_MAX,DELTA_L_MIN,DELTA_U,LAMBDA,u,comp_arco_variavel,Nd,zeta,GDL,NELEM,CONEC,COORDNOS,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,NNOSREST,REST,graph1,DELTA_U_ANTERIOR,fcc_i_anterior,LAMBDA_ANTERIOR)
    
    graph1(n+1,4) = DELTA_L;
    
    delta_u_f = Kcc\fcc;         % Sub incremento devido ao fator de carga (f)
    
    if n == 1
        k_o = (delta_u_f'*delta_u_f)\(fcc'*delta_u_f);
    end
    
    k_n = (delta_u_f'*delta_u_f)\(fcc'*delta_u_f);
    CST = k_n/k_o;
    
    DELTA_LAMBDA = (sqrt(delta_u_f'*delta_u_f))\(DELTA_L);
        
        
    % Sinal do preditor
    pred = DELTA_U'*delta_u_f;
    if pred < 0
        DELTA_LAMBDA = -DELTA_LAMBDA;
    end
    
    DELTA_U_1 = DELTA_LAMBDA*delta_u_f;
    DELTA_U = DELTA_U_1;
    [fcc_i] = forcas_internas_plana(GDL,NELEM,CONEC,COORDNOS,(u+DELTA_U),PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,n,NNOSREST,REST);
    fcc_e = (DELTA_LAMBDA+LAMBDA)*fcc;
    Erro =  fcc_i - fcc_e;                                                  %Vetor das forças internas menos o vetor das forças externas;
    
    
    % INÍCIO DO PROCESSO ITERATIVO
    for i = 1:imax
        
        delta_u_r = -Kcc\Erro;                                              % Sub incremento devido ao fator de carga (r)
        
        delta_lambda = -(DELTA_U_1'*delta_u_f)\(DELTA_U_1'*delta_u_r);
        
        delta_u = delta_u_r + delta_lambda*delta_u_f;                       % Sub incremento final dos deslocamentos
        
        DELTA_U = DELTA_U + delta_u;
        
        DELTA_LAMBDA = DELTA_LAMBDA + delta_lambda;
        
        [fcc_i] = forcas_internas_plana(GDL,NELEM,CONEC,COORDNOS,(u+DELTA_U),PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,n,NNOSREST,REST);   
        
        Erro = fcc_i - (LAMBDA + DELTA_LAMBDA)*fcc;                        %Vetor das forças internas menos o vetor das forças externas;
        
        ACO_F = sqrt(Erro'*Erro);       
        
        if ACO_F <= TOL*sqrt((LAMBDA^2)*(fcc'*fcc))
            break;
        end
        
    end % FIM DO PROCESSO ITERATIVO
    
    if comp_arco_variavel == 1
        DELTA_L = DELTA_L*(Nd/i)^(zeta); 
        if DELTA_L > DELTA_L_MAX 
            DELTA_L = DELTA_L_MAX;
        elseif DELTA_L < DELTA_L_MIN
            DELTA_L = DELTA_L_MIN;
        end
    end
    
    if n > 1 && ACO_F > TOL*sqrt(((LAMBDA)^2)*(fcc'*fcc))
        fprintf('\nNão convergiu no incremento %i',n);
        fprintf('\nLAMBDA = %.2f',LAMBDA);
        fprintf('\nResíduo = %.4f ', ACO_F);
        
        
        % PROCESSO CORRETOR
        DELTA_L = DELTA_L_MIN;
        graph1(n+1,4) = DELTA_L;
        
        delta_u_f = Kcc\fcc;         % Sub incremento devido ao fator de carga (f)
    
        if n == 1
            k_o = (delta_u_f'*delta_u_f)\(fcc'*delta_u_f);
        end

        k_n = (delta_u_f'*delta_u_f)\(fcc'*delta_u_f);
        CST = k_n/k_o;

        DELTA_LAMBDA = (sqrt(delta_u_f'*delta_u_f))\(DELTA_L);

        % Sinal do preditor
        pred = DELTA_U_ANTERIOR'*delta_u_f;
        if pred < 0
            DELTA_LAMBDA = -DELTA_LAMBDA;
        end

        DELTA_U_1 = DELTA_LAMBDA*delta_u_f;
        DELTA_U = DELTA_U_1;
        fcc_i = fcc_i_anterior + Kcc*DELTA_U;
        fcc_e = (DELTA_LAMBDA+LAMBDA_ANTERIOR)*fcc;
        Erro =  fcc_i - fcc_e;  
   
        % PROCESSO ITERATIVO CORRETOR
        for i = 1:imax
       
            delta_u_r = -Kcc\Erro;                                              % Sub incremento devido ao fator de carga (r)

            delta_lambda = -(DELTA_U_1'*delta_u_f)\(DELTA_U_1'*delta_u_r);

            delta_u = delta_u_r + delta_lambda*delta_u_f;                       % Sub incremento final dos deslocamentos

            DELTA_U = DELTA_U + delta_u;

            DELTA_LAMBDA = DELTA_LAMBDA + delta_lambda;
            
            fcc_i = fcc_i + Kcc*delta_u;

            Erro = fcc_i - (LAMBDA + DELTA_LAMBDA)*fcc;                        %Vetor das forças internas menos o vetor das forças externas;

            ACO_F = sqrt(Erro'*Erro);       

            if ACO_F <= TOL*sqrt((LAMBDA^2)*(fcc'*fcc))
                fprintf('\nErro no incremento n = %i -> Houve tentativa de correção',n);
                fprintf('\nTentar com outros parâmetros de comprimento de arco');
                fprintf('\nNão é possível afirmar se houve convergência\n');
                
                figure(2);
                scatter(-graph1(n,3),graph1(n,2),'diamond','k');

                break;
            end

        end % FIM DO PROCESSO ITERATIVO CORRETOR
        
        if i == imax && ACO_F > TOL*sqrt(((LAMBDA)^2)*(fcc'*fcc))
            fprintf('\nNão foi possível a correção no incremento %i',n);
            fprintf('\nTentar com outros parâmetros de comprimento de arco');
            fprintf('\nResíduo = %.4f \n', ACO_F);
        end
        
        % FIM DO PROCESSO CORRETOR
        
    end
    
    u = u + DELTA_U;
    LAMBDA = LAMBDA+DELTA_LAMBDA;
    
    DELTA_U_ANTERIOR = DELTA_U;
    fcc_i_anterior = fcc_i;
    LAMBDA_ANTERIOR = LAMBDA;
    
    
end