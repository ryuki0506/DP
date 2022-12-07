add_library(greetings SHARED field.cpp culc.cpp output.cpp)

target_include_directories(greetings PUBLIC ${PROJECT_SOURCE_DIR}/include)