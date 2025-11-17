cmake -B build
cmake --build build

# intersection test
# ./build/tests/intersection_tests_exe

# bvh test
# ./build/tests/bvh_tests_exe

# refraction
# ./build/src/renderer data/cbox_no_light_refract.json -o cbox_no_light_refract.exr

# reflection
# ./build/src/renderer data/cbox_no_light.json -o cbox_no_light.exr

# area light
./build/src/renderer data/cbox.json -o cbox.exr

# environment light
./build/src/renderer data/sphere_direct.json -o sphere_direct.exr