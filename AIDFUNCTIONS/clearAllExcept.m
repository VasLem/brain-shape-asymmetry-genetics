function test = clearAllExcept(test,in)
    if isempty(test), return; end
    ind = [];
    for i=1:length(in)
        ind = [ind find(strcmp(test,in{i}))];
    end
    test = test(setdiff(1:length(test),ind));
    %for i=1:1:length(test),eval(['clear ' test{i}]);end
end