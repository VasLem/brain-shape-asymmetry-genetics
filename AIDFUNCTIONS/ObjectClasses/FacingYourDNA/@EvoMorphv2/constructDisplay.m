function constructDisplay(obj)
               obj.Fig = figure('Tag','EvoMorph v3.0 copyright Peter Claes',...
                                'UserData',obj, ...
                                'Color',[1 1 1], ...
                                'CloseRequestFcn',@figure_close_Callback);
               set(obj.Fig,'Position',[ 8   436   676   734]);
               set(obj.Fig,'Menubar','none');
               set(obj.Fig,'WindowKeyPressFcn',@KeyPress);
               obj.sError = subplot(4,2,[1 2]);
               set(obj.sError,'xlim',[0 obj.MaxGenerations],'ylim',[0 obj.MeanScore]);
               xlabel(obj.sError,'Generation')
               ylabel(obj.sError,'Fitness Value');
               %obj.sError.FontSize = 20;
               hold on;
               obj.sDiversity = subplot(4,2,[3 4]);
               set(obj.sDiversity,'xlim',[0 obj.MaxGenerations],'ylim',[0 obj.Diversity*2]);
               xlabel(obj.sDiversity,'Generation')
               ylabel(obj.sDiversity,'Diversity');
               %obj.sDiversity.FontSize = 20;
               hold on;
               obj.sScores = subplot(4,2,[5 7]);
               xlabel(obj.sScores,'Score')
               ylabel(obj.sScores,'Nr Individuals');
               set(obj.sScores,'xlim',[obj.BestScore obj.WorstScore],'ylim',[0 obj.PopSize/2]);
               title(obj.sScores,'Score Histogram');
               %obj.sScores.FontSize = 20;
               hold off;
               obj.sScaling = subplot(4,2,[6 8]); 
               xlabel(obj.sScaling,'Raw Scores');
               set(obj.sScaling,'xlim',[0 obj.MeanScore]);
               ylabel(obj.sScaling,'Expectation');
               title(obj.sScaling,'Scaling');
               %obj.sScaling.FontSize = 20;
               hold off;
               obj.vFace = viewerLite;
               obj.ResFace.Axes = obj.vFace.RenderAxes;
               obj.ResFace.Visible = true;obj.ResFace.Selected = true;
               obj.vFace.SceneLightVisible = true;obj.vFace.SceneLightLinked = true;
               set(obj.vFace.Figure,'Position',[140.0000   16.1538  148.6000   43.8462]);
               obj.ResFace.SingleColor = [0.8 0.8 0.8];
               obj.ResFace.Material = 'Dull';
               %set(obj.vFace.Toolbar.light_toggle,'State','on');
               %set(obj.vFace.Toolbar.link_toggle,'State','on');
               obj.PopLM = LMObj('Vertices',obj.Pop(:,obj.vAxes)','Value',obj.Scores);
               obj.PopLM.ColorMode = 'Indexed';
               obj.vPop = viewerLite;
               set(obj.vPop.Figure,'Position',[140.0000   64.9231   80.2000   22.8462]);
               obj.vPop.AxesVisible = true;
               obj.vPop.AxesGrid = true;
               obj.vPop.AxesBox = true;
               obj.vPop.AxesXColor = [0 0 0];
               obj.vPop.AxesYColor = [0 0 0];
               obj.vPop.AxesZColor = [0 0 0];
               obj.PopLM.Axes = obj.vPop.RenderAxes;
               obj.PopLM.Visible = true;obj.PopLM.Selected = true;
               set(obj.vPop.RenderAxes,'clim',[0 obj.WorstScore]);
               set(obj.vPop.RenderAxes,'xlim',[obj.GlobalRange(obj.vAxes(1),1) obj.GlobalRange(obj.vAxes(1),2)]);
               set(obj.vPop.RenderAxes,'ylim',[obj.GlobalRange(obj.vAxes(2),1) obj.GlobalRange(obj.vAxes(2),2)]);
               set(obj.vPop.RenderAxes,'zlim',[obj.GlobalRange(obj.vAxes(3),1) obj.GlobalRange(obj.vAxes(3),2)]);
               %set(obj.vPop.Toolbar.colorbar,'State','on');
               colormap(obj.vPop.RenderAxes,'jet');
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