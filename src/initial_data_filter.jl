# This is the code that I used to get the smaller tif file from the much larger source file. Accessed 2023-02-06.
# You can download the source file here: https://www.sciencebase.gov/catalog/item/6275fd21d34e8d45aa6e2415
# Direct download link: https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/historical/n48w123/USGS_13_n48w123_20220505.tif
@time begin
    using GeoArrays
end
ga = GeoArrays.read("data/USGS_13_n48w123_20220505.tif")
federalway_southeast_corner_index = indices(ga, [-122.2714, 47.25107])
federalway_northwest_corner_index = indices(ga, [-122.4183, 47.36334])
sub_ga = ga[federalway_northwest_corner_index[1]:federalway_southeast_corner_index[1], federalway_northwest_corner_index[2]:federalway_southeast_corner_index[2]]
GeoArrays.write("data/federal_way_elevation.tif", sub_ga)
