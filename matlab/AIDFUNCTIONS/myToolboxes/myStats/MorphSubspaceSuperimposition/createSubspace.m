function QF = createSubspace(F)
         % Trivial column scaling first, if ORTH.m is used later
         for i=1:size(F,2),
             normi=norm(F(:,i),inf);
             %Adjustment makes tol consistent with experimental results
             if normi > eps^.981
               F(:,i)=F(:,i)/normi;
               % Else orth will take care of this
             end
         end
         QF = orth(F);
end