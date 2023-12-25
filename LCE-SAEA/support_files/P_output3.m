function P_output3(Population,Algorithm,Problem,M,Run,p,IGDvalue,FEvalue)

p = p*100;
num2str(p)
d = size(Population,2);
FunctionValue = P_objective1('value',Problem,M,Population);

if(strcmp(Algorithm, 'cRVEA'))
    FunctionValue = FunctionValue(:,1:end - 1);
end
TrueValue = P_objective1('true',Problem,M,1000);

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue = FunctionValue(NonDominated,:);
GDvalue = GD(FunctionValue,TrueValue);
IGDend = IGD(FunctionValue,TrueValue);
PF=size(FunctionValue,1);

if(M == 2)
    figure;
    Plot2D(TrueValue, FunctionValue, 'ro');
%     hold on;
%     scatter(FunctionValue(:,1),FunctionValue(:,2),"green");
%     hold on;
%     scatter(FunctionValue1(:,1),FunctionValue1(:,2),"yellow");
%     hold off;
    legend('All solutions','Nondominated solutions');
end



eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(p),'_',num2str(M),'_',num2str(Run),' Population FunctionValue  GDvalue IGDend IGDvalue FEvalue',])
end


