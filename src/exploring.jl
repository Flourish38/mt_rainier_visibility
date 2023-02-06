@time begin
    @time using GeoArrays
    @time using Plots
    @time federalway_elevation = GeoArrays.read("data/federal_way_elevation.tif")
    @time plot(federalway_elevation)
end

