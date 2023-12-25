% A local correlation estimation surrogate-assisted bi-objective evolutionary
% algorithm for heterogeneous objectives
%------------------------------- Reference --------------------------------
% Chenyan Gu, Handing Wang*, A Local Correlation Estimation Surrogate-
% Assisted Bi-Objective Evolutionary Algorithm for Heterogeneous Objectives,
% Applied Soft Computing, vol.151, pp.111175, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.

% This code is written by Chenyan Gu.
% Email: gugugcy@163.com




Problems={'UF2'};
ps = {0.5};%,1,1.5,2,2.5,3
addpath support_files;
addpath Metric;

for Prob = 1:length(Problems)
    
    Problem=Problems{Prob};
     
    for p = 1:length(ps)
        sheld = ps{p};
        for Run = 1:20
            close all
            LCE('ver',Problem,Run,sheld);
            % If aff, change run_K_RVEA.m please.
        end
        
    end

end
delete(gcp("nocreate"));

function LCE(Algorithm,Problem,Run,p)
    warning off
    
    %% Input: 
    M = 2; % number of objectives
    [~,D] = P_settings('RVEA',Problem,M);
    %number of variables
    no_var = D;
    
    Bounds = [ones(1,no_var);zeros(1,no_var)]; % number of variables and their bounds

    latency = 5; % Latency value
    Max_FE_ex = 200; 
    Max_FE_nex = Max_FE_ex*latency;

                                                                                                             
    THETA = 5.*ones(M+1,D);
    empty_ref = 0;
    id_ex = 2;
    id_nex = 1;

    FE_ex = 0; FE_nex = 0; 
    A = [];
    %A_nex = []; A_ex = [];
    Nex = [];  % Ach
    Ex = [];   % Aex
    Ass = [];  % Aco
    AA = [];

    iter =1; 
    %% Step-1: Generate the initial data
    if no_var>11
        popsize = 110;
    else
        popsize = 10*no_var;
    end
    P = generate_initial_data(Bounds,popsize);

    P = unique(P,'rows');
 
    F_nex = evaluate_least_expensive_obj(P,Problem,id_nex,D);
    F_ex = evaluate_most_expensive_obj(P,Problem,id_ex,D);
    
    Dvalue = F_ex-F_nex;
    Ass = [Ass;[P,Dvalue]];

    FE_nex =FE_nex+popsize;
    FE_ex = FE_ex+popsize;
    FEvalue(iter)= popsize;
    IGDvalue(iter) = P_output2(P,Problem,M);
    Nex = [Nex;[P,F_nex]];
    Ex = [Ex;[P,F_ex]];
    A= [select_solutions_for_archive(P,F_ex,F_nex,id_ex,id_nex)];
    AA = [AA;A];
    

    X_nex  = P;

    s = 0;

     %% Step2
    while (FE_nex<Max_FE_nex && FE_ex<=Max_FE_ex-5)

        % build Kriging for ex,ch,co

        [~,index] = unique( roundn(Nex,-4) ,'rows');
        PopDec1 = Nex(index,1:no_var);
        PopObj1 = Nex(index,1+no_var);
        dmodel     = dacefit(PopDec1,PopObj1,'regpoly1','corrgauss',THETA(1,:),1e-5.*ones(1,D),100.*ones(1,D));
        model_nex   = dmodel;
        THETA(1,:) = dmodel.theta;


        [~,index] = unique( roundn(Ex,-4) ,'rows');

        PopDec2 = Ex(index,1:no_var);
        PopObj2 = Ex(index,1+no_var);
        dmodel     = dacefit(PopDec2,PopObj2,'regpoly1','corrgauss',THETA(2,:),1e-5.*ones(1,D),100.*ones(1,D));
        model_ex   = dmodel;
        THETA(2,:) = dmodel.theta;        

        [~,index] = unique( roundn(Ass,-4) ,'rows');

        PopDec3 = Ass(index,1:no_var);
        PopObj3 = Ass(index,1+no_var);
        dmodel  = dacefit(PopDec3,PopObj3,'regpoly1','corrgauss',THETA(3,:),1e-5.*ones(1,D),100.*ones(1,D));
        model_ass = dmodel;
        THETA(3,:) = dmodel.theta; 
     

        [X_ex,~] = run_K_RVEA(model_ex,model_nex,Bounds,AA,id_ex,id_nex,empty_ref,FE_ex,Max_FE_ex);
        s = s+size(X_ex,1);
        
        F_ex = evaluate_most_expensive_obj(X_ex,Problem,id_ex,no_var);
        F_nex = evaluate_least_expensive_obj(X_ex,Problem,id_nex,no_var);
        FE_ex = FE_ex+size(X_ex,1);
        FE_nex = FE_nex+size(X_ex,1);

        Ex = [Ex;[X_ex,F_ex]];
        Ass = [Ass;[X_ex,F_ex-F_nex]];
        A = [select_solutions_for_archive(X_ex,F_ex,F_nex,id_ex,id_nex)];
        AA = [AA;A];

        PopDec = AA(:,1:no_var);
        PopObj = AA(:,1+no_var:end);
        %% Select promising solutions
        Nondominate = P_sort(PopObj,'half');
        First = find(Nondominate~= Inf);
        PopObj = PopObj(First,:);
        PopDec = PopDec(First,:);
        [lia,locb] = ismember(X_ex,PopDec,'rows');
        index_r = find(lia==0);
        X_add = X_ex(index_r,:);
        index = find(lia~=0);
        X_ex = X_ex(index,:);

        F_nex = F_nex(index,:);
        F_ex = F_ex(index,:);

        N1=(latency-1);

        uniformity = ones(1,size(X_ex,1));
        cub = p*ones(1,no_var);

        %% Local soo
        for i = 1:size(X_ex,1) 
            ub = min((X_ex(i,:)+cub(1,:)),Bounds(1,:));
            lb = max(Bounds(2,:),(X_ex(i,:)-cub(1,:)));
            X = (repmat(ub-lb,N1,1).*lhsamp(N1,D)+repmat(lb,N1,1));
            PopNew1 = X;
            V_nex = evaluate_least_expensive_obj(PopNew1,Problem,id_nex,no_var);
            V_ex = modelvalue(PopNew1,V_nex,model_ass,model_ex);
            FE_nex = FE_nex + size(PopNew1,1);
            V_nex = [V_nex;F_nex(i,:)];
            PopNew1  = [PopNew1;X_ex(i,:)];
            A = [PopNew1,V_nex,[V_ex;F_ex(i,:)]];
            %calculate U
            [a,b,uniformity(i)] = judge(A,no_var,1);
            
            mut_strength = 0.1;
            
            parent = PopNew1;

            record = [];
            child = [];
            un = uniformity(i);

            msoo = latency;

            for n = 1:msoo

                Bound = [ub;lb];
                if b == 0 || un>1
                    [parent,mut_strength,V_nex,Nex] = kill_bad_mnl(parent,V_nex,mut_strength,Problem,id_nex,Bound,'descend',Nex);

                    FE_nex = FE_nex+size(parent,1);
                else

                    [parent,mut_strength,V_nex,Nex] = kill_bad_mnl(parent,V_nex,mut_strength,Problem,id_nex,Bound,'ascend',Nex);
                    FE_nex = FE_nex+size(parent,1);

                    ub = min(ub+cub,Bounds(1,:));
                    lb = max(lb-cub,Bounds(2,:));
                end
                
               
                F_s = modelvalue(parent,V_nex,model_ass,model_ex);

                
                record = [record;[V_nex,F_s]];
                child = [child;parent];
                if std(V_nex,1,1)<1e-5
                    break;
                end
            end

            [~,index] = unique(child,"rows");
            child = child(index,:);
            record = record(index,:); 

            Func = record;
            index = P_sort(Func,'first')==1;
            popnew = child(index,:);

            if size(popnew,1)>=5
                index = randperm(size(popnew,1),5);
                parent = popnew(index,:);
            else
                parent = popnew;
            end

            [V_ex,Ex,FE_ex] = evaluate(Ex,parent,Problem,no_var,FE_ex);
            V_nex = evaluate_least_expensive_obj(parent,Problem,id_nex,no_var);
               
            d_nex = [F_nex(i,:);V_nex];
            d_ex = [F_ex(i,:);V_ex];


            A = [select_solutions_for_archive(parent,V_ex,V_nex,id_ex,id_nex)];
            
            Dvalue = V_ex-V_nex;
            Ass = [Ass;[parent,Dvalue]];
            AA = [AA;A];

        end  

        Pop = Ex(:,1:no_var);
        iter = iter+1;
        IGDvalue(iter)=P_output2(Pop,Problem,M);
        FEvalue(iter) = FE_ex;

    end

    can_A = Ex(:,1:no_var);

    %P_output3(can_A,Algorithm,Problem,M,Run,p,IGDvalue,FEvalue);

end

function X = generate_initial_data(Bounds,sample_size)
    no_var = size(Bounds,2);
    Xn = lhsdesign(sample_size,no_var);
    ub = Bounds(1,:); lb = Bounds(2,:);
    X = bsxfun(@plus,lb,bsxfun(@times,Xn,(ub-lb)));
end


function A = select_solutions_for_archive(P,F_exp,F_nex,id_ex,id_nex)
    S = size(P,1);  
    F = zeros(S,2);
    F(:,id_ex) = F_exp;
    F(:,id_nex) = F_nex;
    A = [P,F];
    
end  


function [f,Ex,FE_ex]=evaluate(Ex,X_ex,Problem,no_var,FE_ex)
    f=ones(size(X_ex,1),1)*1000;
    pop = Ex(:,1:no_var);
    
    for i = 1:size(X_ex,1)
        flag =0;
        x = X_ex(i,:);
        [lia,locb] = ismember(x,pop,'rows');
        if lia==1
            f(i) = Ex(locb,no_var+1);
            flag = 1;
        end
        if flag == 0
            ff=evaluate_most_expensive_obj(x,Problem,2,no_var);
            
            f(i)=ff;
            Ex = [Ex;[x,ff]];
            FE_ex=FE_ex+1;
        end
    end
end