function V_ex = modelvalue(Pop,V_nex,model_ass,model_ex)
%co-surrogate model
        N1 = size(Pop,1);  
        F_c = zeros(N1,1);%co-surrogate mean
        F_s = zeros(N1,1);
        MSEd = zeros(N1,1);

        for j = 1:N1
            x_tr = Pop(j,:);
            F_c(j) = predictor(x_tr, model_ass);
            [F_s(j),~,MSEd(j)] = predictor(x_tr,model_ex);
        end


        uu = F_s+MSEd;
        ll = F_s-MSEd;

        V_ex = F_c+V_nex;
        index1= find(V_ex<ll);
        V_ex(index1,1) = ll(index1,1);

        index2=find(V_ex>uu);
        V_ex(index2,1) = uu(index2,1);
end