function  [] = writeKML_GoogleEarth(Filename, latitude, longitude, altitude, ColorPoints)

%--- header --------------------------------------------------------------%
fod = fopen([Filename '.kml'], 'wb');
fprintf(fod, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fod, '<kml xmlns="http://earth.google.com/kml/2.1">\n');
fprintf(fod, '<Document>\n');
%-------------------------------------------------------------------------%

% create my_style1
fprintf(fod, '<Style id="my_style1">\n');
fprintf(fod, '\t<IconStyle>\n');
fprintf(fod, ['\t\t','<color>ffffffff</color>','\n']);
% fprintf(fod, '\t\t<colorMode>random</colorMode>'); % normal or random
fprintf(fod, '\t\t<scale>0.5</scale>\n');
fprintf(fod, '\t\t<Icon>\n');
fprintf(fod, '\t\t\t<href>http://maps.google.com/mapfiles/kml/pal4/icon26.png</href>\n');
fprintf(fod, '\t\t</Icon>\n');
fprintf(fod, '\t</IconStyle>\n');
fprintf(fod, '</Style>\n\n');
%-------------------------------------------------------------------------%

%--- Write the information for each point --------------------------------%
for n=1:length(altitude)
    fprintf(fod, '<Placemark>\n');
    %     fprintf(fod, '<name>%d</name>\n',n);
    %     fprintf(fod, '<name>%d</name>\n',n-6); % in case of test6_repl.ubx file
    
    fprintf(fod, '\t<styleUrl>#my_style1</styleUrl>\n');
    
    fprintf(fod, '\t<Point>\n');
    fprintf(fod, '\t<text>ccc</text>\n');
    fprintf(fod, '\t\t<coordinates>\n');
    fprintf(fod,'\t\t\t%.20f, %.20f, %.20f\n ',longitude(n),latitude(n),altitude(n));
    fprintf(fod, '\t\t</coordinates>\n');
    fprintf(fod, '\t</Point>\n');
    fprintf(fod, '</Placemark>\n');
end
%-------------------------------------------------------------------------%


fprintf(fod, '\n\n</Document>\n');
fprintf(fod, '</kml>\n');
fclose(fod);
