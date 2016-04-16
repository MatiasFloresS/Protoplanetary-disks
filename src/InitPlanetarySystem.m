function cant_planets = InitPlanetarySystem(filename)
    global system;
    sprintf(filename);
    cant_planets = FindNumberOfPlanets(filename);
    sprintf ('%d planet(s) found.\n', cant_planets);
    AllocPlanetsSystem(cant_planets);
    system{1} = cant_planets;
    information_planets(filename,cant_planets);
    return
end