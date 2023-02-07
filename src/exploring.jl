@time begin
    @time using GeoArrays
    @time using Plots
    coords = GeoArrays.coords
    using PlyIO
    using Statistics
    using LinearAlgebra
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

points_matrix = reduce(vcat, transpose.(points))

mean_point = vec(mean(points_matrix, dims=1))
mean_point_norm = mean_point ./ sqrt(sum(mean_point .^ 2))
target = [0, 0, 1.]

function rotation_matrix(from, to)
    temp = from + to
    return 2(temp*temp')/(temp'*temp) - I
end

R = rotation_matrix(mean_point_norm, target)

transformed_points = points_matrix * R

ply = Ply()
push!(ply, PlyElement("vertex", 
    ArrayProperty("x", transformed_points[:, 1]),
    ArrayProperty("y", transformed_points[:, 2]),
    ArrayProperty("z", transformed_points[:, 3])
))

save_ply(ply, "data/federal_way_true.ply")