function screenPoint2currentPoint(obj)
         scrn_pt = get(0, 'PointerLocation');
         set(obj.Figure,'units','pixels')
         loc = get(obj.Figure, 'Position');
         % We need to compensate for an off-by-one error:
         pt = [scrn_pt(1) - loc(1) + 1, scrn_pt(2) - loc(2) + 1];
         set(obj.Figure,'CurrentPoint',pt);
end