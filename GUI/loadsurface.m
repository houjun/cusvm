function surface = loadsurface(filename)

fp = fopen(filename,'r');
surface.numLat = fscanf(fp,' %d ',1);
surface.numLon = fscanf(fp,' %d ',1);
surface.lati = fscanf(fp,'%f',surface.numLat);
surface.loni = fscanf(fp,'%f',surface.numLon);
surface.raster = fscanf(fp,'%f',[surface.numLon,surface.numLat])';

fclose(fp);
