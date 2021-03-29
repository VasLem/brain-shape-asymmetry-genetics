
function f = showSignificanceEvolution(brainSurface, nSamplesPerPick, perm, title)
f = uifigure('HandleVisibility', 'on', 'Position', [60 60 350 275],'Name',title);
ax = uiaxes(f, 'Position',[60 120 240 240]);
 axis(ax, 'equal');
showing= renderBrainSurface(brainSurface,perm(1,:));
showing.RenderAxes = ax;
sld = uislider(f,'ValueChangedFcn',@(sld,event) updateVertexValue(sld,showing, perm));
bt = uibutton(f,'state', 'Text', 'Play', 'Value', false, 'Position', [40, 70, 50,40]);
t = timer('Period',0.2,'ExecutionMode', 'fixedRate');
t.TimerFcn =  @(ax, event) PlayOrStop(t, f, bt, sld, showing, perm);

sld.MajorTickLabels = string(nSamplesPerPick);
sld.MajorTicks = 1:size(nSamplesPerPick,2);
sld.Limits =  [1 length(nSamplesPerPick)];
start(t);

function updateVertexValue(sld, showing, perm)
sld.Value = round(sld.Value);
showing.VertexValue = perm(sld.Value,:);
end


function PlayOrStop(timer, fig, bt, sld, showing, perm)
if ~isvalid(fig)
    stop(timer);
end
if bt.Value
    sld.Value = max(1, mod(sld.Value+ 1, size(perm,1)+1));
    updateVertexValue(sld,showing, perm);
end
end
end