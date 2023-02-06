@time begin
    @time using GeoArrays
    @time using Plots
    coords = GeoArrays.coords
    using Geodesy
    @time federalway_elevation = GeoArrays.read("data/federal_way_elevation.tif")
    @time plot(federalway_elevation)
end

for count in 1:100 
    ll = rand(coords(federalway_elevation))
    ll_indices = indices(federalway_elevation, ll, Center())
    a = federalway_elevation[ll_indices..., 1]
    a2 = federalway_elevation[ll_indices]
    if a != a2[1]
        print(a, "\t", a2, "\t", count)
        break
    end
end

begin
    ll = rand(coords(federalway_elevation))
    ll_indices = indices(federalway_elevation, ll, Center())
    a = federalway_elevation[ll_indices..., 1]
    x_lla = LLA(ll..., a)
    x_ecef = ECEF(x_lla, wgs84)
end

@time points = vec([ECEFfromLLA(wgs84)(LLA(ll..., federalway_elevation[indices(federalway_elevation, ll, Center())..., 1])) for ll in coords(federalway_elevation)])