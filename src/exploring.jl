@time begin
    @time using GeoArrays
    using PlyIO
    using Statistics
    using LinearAlgebra
    using Geodesy
    federal_way_elevation = GeoArrays.read("data/federal_way_elevation.tif")
    mt_rainier_elevation = GeoArrays.read("data/mt_rainier_elevation.tif")
end

for count in 1:100 
    ll = rand(coords(federal_way_elevation))
    ll_indices = indices(federal_way_elevation, ll, Center())
    a = federal_way_elevation[ll_indices..., 1]
    a2 = federal_way_elevation[ll_indices]
    if a != a2[1]
        print(a, "\t", a2, "\t", count)
        break
    end
end

begin
    ll = rand(coords(federal_way_elevation))
    ll_indices = indices(federal_way_elevation, ll, Center())
    a = federal_way_elevation[ll_indices..., 1]
    x_lla = LLA(ll..., a)
    x_ecef = ECEF(x_lla, wgs84)
end

@time points = vec([ECEFfromLLA(wgs84)(LLA(ll..., federal_way_elevation[indices(federal_way_elevation, ll, Center())..., 1])) for ll in coords(federal_way_elevation)])

points_matrix = reduce(vcat, transpose.(points))

ll_center = coords(federal_way_elevation, size(federal_way_elevation)[1:2] .÷ 2)
center_surface_normal = [(ECEFfromLLA(wgs84)(LLA(ll_center..., 1)) - ECEFfromLLA(wgs84)(LLA(ll_center..., 0)))...]

#=
# This does not actually compute the surface normal... obviously.
mean_point = vec(mean(points_matrix, dims=1))
mean_point_norm = mean_point ./ sqrt(sum(mean_point .^ 2))
=#

function rotation_matrix(from, to)
    temp = from + to
    return 2(temp*temp')/(temp'*temp) - I
end

R1 = rotation_matrix(center_surface_normal, [0, 0, 1.])

transformed_points = points_matrix * R1
offsets_1 = [mean(transformed_points[:, 1:2], dims=1)..., minimum(transformed_points[:, 3])]
for i in 1:3
    transformed_points[:, i] .-= offsets_1[i]
end

center_point_index = size(federal_way_elevation)[1]÷2 + size(federal_way_elevation)[1]*(size(federal_way_elevation)[2]÷2)
center_point_flat = vcat(transformed_points[center_point_index, 1:2], 0)
r1_east_vector = vcat(transformed_points[center_point_index + 1, 1:2], 0) .- center_point_flat
r1_east_vector_norm = r1_east_vector ./ norm(r1_east_vector)
#=
# This does not calculate the direction of north... obviously.
center_north_point = vcat(transformed_points[size(federal_way_elevation)[1] ÷ 2, 1:2], 0)
center_north_point_norm = center_north_point ./ norm(center_north_point)
=#
R2 = rotation_matrix(center_north_point_norm, [0, -1., 0]) * [1 0 0; 0 -1 0; 0 0 -1]

transformed_points = transformed_points * R2

offsets_2 = [mean.(extrema(transformed_points[:, 1:2], dims=1))..., 0]

for i in 1:3
    transformed_points[:, i] .-= offsets_2[i]
end


ply = Ply()
push!(ply, PlyElement("vertex", 
    ArrayProperty("x", transformed_points[:, 1]),
    ArrayProperty("y", transformed_points[:, 2]),
    ArrayProperty("z", transformed_points[:, 3])
))

save_ply(ply, "data/federal_way_true.ply")