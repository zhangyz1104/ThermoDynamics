project(ThermoDynamics)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES lib/main.c lib/global_defs.h lib/Instructions.c lib/prototypes.h lib/Geodynamics.c lib/Crystallization.c)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} m)