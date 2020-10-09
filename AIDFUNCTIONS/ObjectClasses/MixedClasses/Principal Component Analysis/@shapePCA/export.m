function export(obj)
        path = uigetdir(pwd,'Export Folder');
        if path==0, return; end
        cd(path);
        exportWavefront(obj.Average,'Model');
        % Exporting EigenVectors
        EigVec = obj.EigVec;
        save('Model.EigVec','-ASCII','EigVec');
        % Exporting EigenValues
        EigVal = obj.EigVal;
        save('Model.EigVal','-ASCII','EigVal');
        % Exportin TCoeff
        Tcoeff = obj.Tcoeff;
        save('Model.Tcoeff','-ASCII','Tcoeff');
end