# This is the code that I used to get the smaller tif file from the much larger source file. Accessed 2023-02-06.
# You can download the source file here: https://www.sciencebase.gov/catalog/item/6275fd21d34e8d45aa6e2415
# Direct download link: https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/historical/n48w123/USGS_13_n48w123_20220505.tif
@time begin
    using GeoArrays
end

begin
    ga = GeoArrays.read("data/USGS_13_n48w123_20220505.tif")
    federalway_southeast_corner_index = indices(ga, [-122.2714, 47.25107])
    federalway_northwest_corner_index = indices(ga, [-122.4183, 47.36334])
    sub_ga = ga[federalway_northwest_corner_index[1]:federalway_southeast_corner_index[1], federalway_northwest_corner_index[2]:federalway_southeast_corner_index[2]]
    GeoArrays.write("data/federal_way_elevation.tif", sub_ga)
end

begin
    @time ga = GeoArrays.read("data/USGS_13_n47w122_20220919.tif")
    mt_rainier_index = indices(ga, collect(coords(ga, Tuple(argmax(ga))[1:2])))
    checked_indices = Set([mt_rainier_index])
    indices_to_check = Set([mt_rainier_index])
    max_x = 0
    min_x = size(ga)[1]+1
    max_y = 0
    min_y = size(ga)[2]+1
    n = 0
    start_time = time()
    while !isempty(indices_to_check)
        n += 1
        if (n % 100000 == 0)
            println(time()-start_time, "\t", (min_x, min_y), "\t", (max_x, max_y), "\t", n)
        end
        index = pop!(indices_to_check)
        push!(checked_indices, index)
        index_elevation = ga[index..., 1]
        # Makes sure the elevation is high enough (~ tree line) and that the hill is sloping down towards federal way (otherwise, it's occluded!)
        if index_elevation >= 2000 && (ga[index[1]-1, index[2], 1] + ga[index[1], index[2]-1, 1] < 2 * index_elevation)
            max_x = max(index[1], max_x)
            min_x = min(index[1], min_x)
            max_y = max(index[2], max_y)
            min_y = min(index[2], min_y)

            for shift in [[1, 0], [-1, 0], [0, 1], [0, -1]]
                new_index = index .+ shift
                if !(new_index in checked_indices)
                    push!(indices_to_check, new_index)
                end
            end
        end
    end
    println(time()-start_time, "\t", (min_x, min_y), "\t", (max_x, max_y), "\t", n)
    sub_ga = ga[min_x:max_x, min_y:max_y]
    GeoArrays.write("data/mt_rainier_elevation.tif", sub_ga)
end
