function fig = figure_by_name(nm)
% Returns the handle of a figure object when passed the name
% Usage: fig28a = figure_by_name('Figure 2.8(a)')
[flag, fig] = figflag(nm);
if flag == 0
    fig = figure;
end
set(fig, 'NumberTitle', 'off', 'Name', nm);
figure(fig);

