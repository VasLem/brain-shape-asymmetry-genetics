function out = getComputerID
         name = getComputerName;
         switch name
             case 'morpht7600' % WIN PETER (16)
                 out = 0;
             case 'micdp02'% PETER UBUNTU (36)
                 out = 2;    
             case 'micpro.uz.kuleuven.ac.be' % MAC JAIRUI (12)
                 out = 2;
             case 'micpro02.uz.kuleuven.ac.be' % MAC DOMINIQUE (12)
                 out = 3;
             case 'micpro03.uz.kuleuven.ac.be' % MAC PHILIP (12)
                 out = 4;
             case 'micapp.uz.kuleuven.ac.be' % MAC INE (12)
                 out = 5;
             case 'radmevislab.uz.kuleuven.ac.be' % MAC BART (6)
                 out = 6;
             case 'micim07.uz.kuleuven.ac.be' % IMAC (4)
                 out = 7;
             case 'naphat.uz.kuleuven.ac.be' % OLD DOMINIQUE (12)
                 out = 8;
             case 'mircmevislab.uz.kuleuven.ac.be' % STUDENT COMPUTER (8)
                 out = 9;
             case 'micim02.uz.kuleuven.ac.be' % STUDENT COMPUTER (4)
                 out = 10;
             case 'micim10.uz.kuleuven.ac.be' % OLD JAAP (4)
                 out = 11;
             case 'micim11.uz.kuleuven.ac.be' % DOROTHY (4)
                 out = 12;
             case 'micim12.uz.kuleuven.ac.be' % DZEMILA (4)
                 out = 13;
             case 'micim08.uz.kuleuven.ac.be' % DAVID (4)
                 out = 14;    
             case 'micpro01.uz.kuleuven.ac.be' % MAC PETER (12)
                 out = 15;
             case 'micim13.uz.kuleuven.ac.be' % SOFIE (4)
                 out = 16;
             case 'micbb01' % CLUSTER1 (24)
                 out = 17;
             case 'micsd01' % CLUSTER2 (24)
                 out = 18;
             case 'micsm01'
                 out = 20;
             case 'michpz01'
                 out = 21;
             case 'michpz03'
                 out = 22;    
             otherwise
                 out = 19;
         end
end