function r = myGetRank(Data)
         data_sorted = sort(Data);
         [~, rnk] = ismember(Data,data_sorted);
         r = rnk(:);
end