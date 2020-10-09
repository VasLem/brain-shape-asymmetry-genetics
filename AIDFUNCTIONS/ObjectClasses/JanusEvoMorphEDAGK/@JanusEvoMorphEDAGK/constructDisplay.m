function constructDisplay(obj)
               obj.Fig = figure('Tag','JanusEvoMorph v1.0 copyright Peter Claes',...
                                'UserData',obj, ...
                                'CloseRequestFcn',@figure_close_Callback);
               set(obj.Fig,'Position',[ 8   436   676   734]);
               set(obj.Fig,'Menubar','none');
               set(obj.Fig,'WindowKeyPressFcn',@KeyPress);
               obj.sError = subplot(4,2,[1 2]);
               set(obj.sError,'xlim',[0 obj.MaxGenerations],'ylim',[0 obj.WorstScore]);
               xlabel(obj.sError,'Generation')
               ylabel(obj.sError,'Fitness Value');
               hold on;
               obj.sDiversity = subplot(4,2,[3 4]);
               set(obj.sDiversity,'xlim',[0 obj.MaxGenerations],'ylim',[0 obj.Diversity*2]);
               %xlabel(obj.sDiversity,'Generation')
               ylabel(obj.sDiversity,'Diversity');
               hold on;
               obj.sScores = subplot(4,2,[5 7]);
               xlabel(obj.sScores,'Score')
               ylabel(obj.sScores,'Nr Individuals');
               title(obj.sScores,'Score Histogram');
               hold off;
               obj.sScaling = subplot(4,2,[6 8]); 
               xlabel(obj.sScaling,'Raw Scores')
               ylabel(obj.sScaling,'Expectation');
               title(obj.sScaling,'Scaling');
               hold off;
               obj.vFace = viewer(obj.ResFace);
               set(obj.vFace.Figure,'Position',[140.2000   61.3077   80.2000   26.5385]);
               obj.ResFace.SingleColor = [0.8 0.8 0.8];
               obj.ResFace.Material = 'Dull';
               set(obj.vFace.Toolbar.light_toggle,'State','on');
               set(obj.vFace.Toolbar.link_toggle,'State','on');
               obj.PopLM = LMObj('Vertices',obj.Population(:,obj.vAxes)','Value',obj.Scores);
               obj.PopLM.ColorMode = 'Indexed';
               obj.vPop = viewer(obj.PopLM);
               set(obj.vPop.Figure,'Position',[140.2000   33.4615   80.2000   22.7692]);
               obj.vPop.AxesVisible = true;
               obj.vPop.AxesGrid = true;
               obj.vPop.AxesBox = true;
               set(obj.vPop.RenderAxes,'clim',[0 obj.WorstScore]);
               set(obj.vPop.RenderAxes,'xlim',[obj.GlobalRange(obj.vAxes(1),1) obj.GlobalRange(obj.vAxes(1),2)]);
               set(obj.vPop.RenderAxes,'ylim',[obj.GlobalRange(obj.vAxes(2),1) obj.GlobalRange(obj.vAxes(2),2)]);
               set(obj.vPop.RenderAxes,'zlim',[obj.GlobalRange(obj.vAxes(3),1) obj.GlobalRange(obj.vAxes(3),2)]);
               set(obj.vPop.Toolbar.colorbar,'State','on');
               % wait untill with pointer over main figure a key is pressed
               while strcmp(obj.State,'Init');
                   disp('press any key to continue');
                   pause(5);
               end
end
function KeyPress(varargin)
        obj = get(varargin{1},'UserData');
        obj.State = 'Run';
end
function figure_close_Callback(varargin)
        obj = get(varargin{1},'UserData');
        obj.State = 'Stop';
        delete(obj.vFace);
        delete(obj.vPop);
        delete(gcf);
end