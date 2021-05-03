classdef spheroid < matlab.mixin.CustomDisplay
%SPHEROID Abstract superclass for oblateSpheroid and referenceSphere
%
%   This class customizes the display for its subclasses, which must define
%   DerivedProperties and DisplayFormat properties.

% Copyright 2012-2015 The MathWorks, Inc.

%-------------------------- Customize disp ----------------------------

properties (Constant, Abstract, Hidden, Access = protected)
    % Cell string listing names "derived" properties, to be displayed with
    % their values suppressed (even for scalar objects).
    DerivedProperties
    
    % Custom format string: '', or a string accepted by the FORMAT
    % function('short', 'longG', ...).  If non-empty, it is used to
    % override the current display formatting for floating point values.
    DisplayFormat
end

properties (Constant, Hidden, Access = private)
    % Use the map:geodesy message catalog to construct the "derived
    % properties header," to enable localization. This message has a
    % constant string, so optimize by caching the message object here and
    % avoiding the need for access each time disp is called on a subclass
    % instance.
    DerivedPropertiesHeaderMessage ...
        = message('map:geodesy:derivedSpheroidPropertiesHeader');
end

methods (Access = protected)
    
    function displayScalarObject(self)
        % (1) Overload the default and suppress display of "derived"
        %     properties, without making them hidden.
        %
        % (2) When non-empty, use the specified DisplayFormat instead of
        %     the current format setting.
        
        header = getHeader(self);
        disp(header);
        
        % Construct a structure that contains the defining properties and
        % their values, and excludes the derived properties.
        definingProperties = setdiff( ...
            properties(self), self.DerivedProperties, 'stable');
        
        values = cellfun(@(propertyName) self.(propertyName), ...
            definingProperties, 'UniformOutput', false);
        
        s = cell2struct(values, definingProperties, 1);
        
        % Manage numerical display format for defining properties.
        if ~isempty(self.DisplayFormat)
            fmt = get(0,'format');
            clean = onCleanup(@() format(fmt));
            format(self.DisplayFormat)
        end
        
        disp(s);
        
        % Display names only for derived properties, taking care to be
        % consistent with the current format spacing.
        looseSpacing = strcmp(get(0,'FormatSpacing'),'loose');
        if ~looseSpacing
            % We want a vertical space before the "derived properties
            % header" in any case.  It comes automatically with 'loose',
            % but needs to be added with 'compact'.
            fprintf('\n')
        end
        fprintf('  %s\n\n', getString(self.DerivedPropertiesHeaderMessage))
        fprintf('    %s\n', self.DerivedProperties{:})
        if looseSpacing
            % Add a trailing newline only when spacing is loose.
            fprintf('\n')
        end
        
        footer = getFooter(self);
        if ~isempty(footer)
            disp(footer);
        end
    end
    
    function header = getHeader(self)
        if ~isscalar(self)
            header = getHeader@matlab.mixin.CustomDisplay(self);
        else
            % Use the map:geodesy message catalog to construct the custom
            % header, to enable localization.
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(self);
            msg = message('map:geodesy:definingSpheroidPropertiesHeader', className);
            header = sprintf('%s\n', getString(msg));
        end
    end
    
end
end
