function myCleanAxes(axh, args)
% based on Sams function adapted by me.
% its meant to clear figure axes to save it for creating nice figures with
% Inkscape/Illustrator
%
% Magdalena Sabat, 05/2024, magdalenajoannasabat@gmail.coom
arguments
    axh {mustBeAxesHandle}
    args.ticklength = 0.02;
    args.tickwidth = 2;
    args.notick = false;
    args.axes = true;
    args.title = true;
end

box off;
if args.notick
    set(axh, 'XTick', [], 'YTick', [], 'ZTick', []);
else
    set(axh, 'TickLength', [args.ticklength, 0], 'linewidth', args.tickwidth);
end
set(axh, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
xlabel(axh,''); ylabel(axh,'');
if ~args.axes
    set(axh, 'visible', 'off');
    set(findall(axh, 'type', 'text'), 'visible', 'on')
end
if ~args.title
    title(axh,'');
end


end
% Custom validation function
function mustBeAxesHandle(axh)
    if ~ strcmp(get(axh,'type'),'axes')
        eid = 'First argument must be axes handle.';
        msg = 'bad argument';
        error(eid,msg)
    end
end