function refMap(lat,lon, name)

ax=worldmap();
setm(ax,'mlabelparallel',-90)
set(ax,'Visible','On')
%land = shaperead('landareas', 'UseGeoCoords', true);
%geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85])
 load coastlines
 plotm(coastlat,coastlon,'k')
hold on
plotm(lat,lon,'ro','MarkerSize',10)
plotm(lat,lon,'r+','MarkerSize',30,'LineWidth',1)
textm(lat,lon,name)
patches = findobj(gca, 'type', 'patch');
set(patches,'Visible','Off')
end