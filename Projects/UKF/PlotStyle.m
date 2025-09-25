classdef PlotStyle
    %PlotStyle  Simple helper to enforce consistent plotting defaults.
    %
    %   style = PlotStyle();                 % use defaults
    %   style = PlotStyle(14,'Arial',1.5);   % custom values
    %
    %   style.apply();   % apply to all new figures/axes
    %
    %   style.apply(fig) % apply to a specific figure (optional)
    %
    %   The properties fs, fn, lw are public so you can read/modify.

    properties
        fs double = 12             % font size
        fn char   = 'Times New Roman'  % font name
        lw double = 2              % default line width
    end

    methods
        function obj = PlotStyle(fs, fn, lw)
            if nargin >= 1 && ~isempty(fs), obj.fs = fs; end
            if nargin >= 2 && ~isempty(fn), obj.fn = fn; end
            if nargin >= 3 && ~isempty(lw), obj.lw = lw; end
        end

        function apply(obj, fig)
            % Apply defaults globally or to a specific figure/axes.
            if nargin < 2 || isempty(fig)
                % global defaults (affects all new figures)
                set(groot,'DefaultAxesFontSize',obj.fs, ...
                          'DefaultAxesFontName',obj.fn, ...
                          'DefaultLineLineWidth',obj.lw, ...
                          'DefaultTextFontSize',obj.fs, ...
                          'DefaultTextFontName',obj.fn);
            else
                % apply to one figure (existing axes)
                ax = findall(fig,'type','axes');
                set(ax,'FontSize',obj.fs,'FontName',obj.fn);
                ln = findall(fig,'type','line');
                set(ln,'LineWidth',obj.lw);
            end
        end
    end
end