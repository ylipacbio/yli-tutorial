project(
    'my',
    ['c', 'cpp'],
    version:'0.1.0',
    default_options: [
        'cpp_std=c++14',
        'c_std=gnu11',
    ],
    meson_version: '>= 0.46.0'
)

my_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
    '-Wno-unused-command-line-argument',
    '-Wno-missing-field-initializers',
    '-Wno-unused-parameter',
    '-Wno-unused-variable',
    '-Wno-unused-const-variable'
    ]
    if (cpp.has_argument(cflag))
        my_flags += cflag
    endif
endforeach

if cpp.get_id() == 'clang'
    foreach cflag: [
        '-Wno-unused-local-typedef',
        '-Wno-c99-extensions']
        if cpp.has_argument(cflag)
            my_flags += cflag
        endif
    endforeach
endif

#if get_option('sse41')
#    if cpp.has_argument('-msse4.1')
#        my_flags += '-msse4.1'
#    else
#        error('Requested SSE4.1, but the C++ compiler does not support -msse4.1')
#    endif
#endif

my_parasail_dep = dependency('parasail-1', required:true, fallback: ['parasail', 'parasail_dep'])

my_lib_deps = [
    my_parasail_dep,
]

my_cpp_sources = files(['my-parasail.cpp'])

my_lib = static_library(
    'my',
    my_cpp_sources,
    install: false,
    dependencies: my_lib_deps,
    cpp_args: my_flags
)

my_main = executable(
    'my', 
    my_cpp_sources,
    install: true,
    dependencies: my_lib_deps,
    include_directories: include_directories('.'),
    link_with: my_lib,
    cpp_args: my_flags
)
