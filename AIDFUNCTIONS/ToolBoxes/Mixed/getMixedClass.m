function C = getMixedClass(diff,Types)
        n = length(Types);
        ClassTypes = zeros(1,n);
        for i=1:1:n
            switch sign(diff(i))
                case 1
                  ClassTypes(i) = 2;
                case -1
                  ClassTypes(i) = 0;
            end
        end
        %C = Types-ClassTypes;
        indexAA = (find(Types==0));
        TAA = length(find(ClassTypes(indexAA)==0));
        CAA = (TAA/length(indexAA))*100;
        indexBB = (find(Types==2));
        TBB = length(find(ClassTypes(indexBB)==2));
        CBB = (TBB/length(indexBB))*100;
        C = ((TBB+TAA)/(length(indexAA)+length(indexBB)))*100;
end