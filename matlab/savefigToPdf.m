function savefigToPdf(fig,path)
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(fig,'RendererMode','manual');
fig.Renderer='Painters';
print(fig,path, '-dpdf');
end

