function(build_hdf5_test)
    get_filename_component(test_name ${ARGV0} NAME_WE)
    add_executable(${test_name} ${ARGV})
    target_link_libraries(${test_name} PUBLIC SZ3 hdf5sz3)
endfunction(build_hdf5_test)

build_hdf5_test(des-w.cpp)
build_hdf5_test(des.cpp)
# build_hdf5_test(src/convertBinToHDF5.cpp)

