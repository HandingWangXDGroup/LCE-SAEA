function [a,b,uniformity] = judge(A,no_var,FE_ex)
    a=0;
    b=0;
    %t1=0.4+0.1*(-0.5*cos(FE_ex*pi/200)+0.5);
    %t2=0.6-0.1*(-0.5*cos(FE_ex*pi/200)+0.5);
    s1 =0;
    s2 = 0;
    f_ex = A(:,no_var+2);
    f_nex = A(:,no_var+1);
    
    [~,index1] = sort(f_nex,'descend');
    r = f_ex(index1);
    for i =1:(length(index1)-1)
        for j = (i+1):length(f_nex)
            if r(i)<r(j) %非支配 
                a=a+1;                                                                                                                                                                                
                s1 = s1+(r(i)-r(j));
            else
                b = b+1;%支配
                s2 = s2+(r(i)-r(j));
            end
        end
    end
    uniformity = a/b;
    %disp(s1);disp(s2);
end