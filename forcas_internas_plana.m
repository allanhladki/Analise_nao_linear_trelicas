% TRELI�A PLANA N�O LINEAR
% AUTOR: ALLAN DJONES COSTA HLADKI
% UNIVERSIDADE TECNOL�GICA FEDERAL DO PARAN�
% GUARAPUAVA - 2021
% Contato: allanhladki@alunos.utfpr.edu.br

% C�LCULO DO VETOR DE FOR�AS INTERNAS
% 

function[fcc_i] = forcas_internas_plana(GDL,NELEM,CONEC,COORDNOS,deslocamentos,PROPELEM,tipo_deformacao,tipo_linear_fisico,SIGMAANT,sigma,normal,n,NNOSREST,REST)

    fcc_i = zeros(GDL,1);
    
    for k = 1: +1 : NELEM 

        % Geometria da barra: comprimento, cossenos diretores
        noi = CONEC(k,1);         %N� inicial do elemento
        nof = CONEC(k,2);         %N� final do elemento

        % Dados iniciais
        dxi = COORDNOS(nof,1) - COORDNOS(noi,1);       % Deltax do elemento (seu comprimento na dire��o x);
        dyi = COORDNOS(nof,2) - COORDNOS(noi,2);       % Deltay do elemento (seu comprimento na dire��o y);
        Le = sqrt(dxi^2 + dyi^2);                      % Comprimento inicial do elemento 

        % Dados ap�s deforma��o uant
        dxf = COORDNOS(nof,1) - COORDNOS(noi,1) + deslocamentos(2*nof-1) - deslocamentos(2*noi-1);       % Deltax do elemento (seu comprimento na dire��o x);
        dyf = COORDNOS(nof,2) - COORDNOS(noi,2) + deslocamentos(2*nof) - deslocamentos(2*noi);           % Deltay do elemento (seu comprimento na dire��o y);
        Lf = sqrt(dxf^2 + dyf^2);                                                                        % Comprimento final do elemento              

        % Pelos cossenos diretores 
        oxf = dxf/Lf;                % Cosseno do �ngulo theta x do vetor;
        oyf = dyf/Lf;                % Cosseno do �ngulo theta y do vetor;

        % Propriedades do elemento
        Ae = PROPELEM(k,1);
        Ee = PROPELEM(k,2);
        poisson = PROPELEM(k,8);
        
        if tipo_deformacao == 1
            eatual = (Lf-Le)/Le;                    % Deforma��o atual do elemento - Engenharia;
        
        elseif tipo_deformacao == 2
            eatual = (Lf^2-Le^2)/(2*Le^2);          % Deforma��o atual do elemento - Green;
        
        elseif tipo_deformacao == 3
            eatual = log(Lf/Le);                    % Deforma��o atual do elemento - Log sem mudan�a de volume;
        
        elseif tipo_deformacao == 4
            eatual = log(Lf/Le);                    % Deforma��o atual do elemento - Log com mudan�a de volume;
            
        end
    
        
        if tipo_linear_fisico == 2
            Eesc = PROPELEM(k,3);                   % M�dulo de elasticidade ap�s o escoamento
            Tesc = PROPELEM(k,4);                   % Tens�o atual de escoamento do elemento
            e_esc = PROPELEM(k,5);                  % Deforma��o de escoamento do elemento
            e_p = PROPELEM(k,6);                    % �ltima deforma��o tal que a tens�o � nula no elemento
            contador = PROPELEM(k,7);               % Contador para c�lculo das tens�es do elemento
            sigmaant = SIGMAANT(k);                 % Tens�o anterior;

            [sigma(k),~,~,~,~,~] = tensao_nlf(Tesc,e_esc,e_p,sigmaant,contador,Ee,Eesc,eatual,n);

            
        elseif tipo_linear_fisico == 1
            [sigma(k)] = tensao_hooke(eatual,Ee);
        end
        
        normal(k) = sigma(k)*Ae;
            
        if tipo_deformacao == 1
            fcc_i_elemento = normal(k)*[-oxf;-oyf;oxf;oyf;];
        elseif tipo_deformacao == 2
            fcc_i_elemento = normal(k)*(Lf/Le)*[-oxf;-oyf;oxf;oyf;];
        elseif tipo_deformacao == 3
            fcc_i_elemento = normal(k)*(Le/Lf)*[-oxf;-oyf;oxf;oyf;];
        elseif tipo_deformacao == 4
            fcc_i_elemento = normal(k)*((Le/Lf)^(2*poisson))*[-oxf;-oyf;oxf;oyf;];
        end      

        % Assemblagem de fcc_i
        MATRIZ_ASSEMBLAGEM = [2*noi-1,2*noi,2*nof-1,2*nof];
      
        for h = 1 : 4
            lin = MATRIZ_ASSEMBLAGEM(h);
            fcc_i(lin) = fcc_i(lin) + fcc_i_elemento(h);
        end
    end     % Fim da cria��o da matriz for�as internas da estrutura
        

    % SISTEMA RESTRINGIDO - CONDI��ES DE CONTORNO   
    for k = 1 : NNOSREST
        no = REST(k,1);
        % Deslocamentos prescritos no vetor fcc
        for m = 2:3
            if REST(k,m) == 1
                GDLR = 2*no-(3-m);
                % Atribui o desloca. pres. na linha GDLR de fcc
                fcc_i(GDLR) = 0;
            end
        end    
    end
    
end
