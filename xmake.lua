set_project("geotiv")
set_version("3.0.0")
set_xmakever("2.7.0")

-- Set C++ standard
set_languages("c++20")

-- Add build options
add_rules("mode.debug", "mode.release")

-- Compiler warnings and flags (matching CMake)
add_cxxflags("-Wall", "-Wextra", "-Wpedantic", "-Wno-reorder")

-- Add global search paths for packages in ~/.local
local home = os.getenv("HOME")
if home then
    add_includedirs(path.join(home, ".local/include"))
    add_linkdirs(path.join(home, ".local/lib"))
end

-- Add devbox/nix paths for system packages
local cmake_prefix = os.getenv("CMAKE_PREFIX_PATH")
if cmake_prefix then
    add_includedirs(path.join(cmake_prefix, "include"))
    add_linkdirs(path.join(cmake_prefix, "lib"))
end

local pkg_config = os.getenv("PKG_CONFIG_PATH")
if pkg_config then
    -- Split PKG_CONFIG_PATH by ':' and process each path
    for _, pkgconfig_path in ipairs(pkg_config:split(':')) do
        if os.isdir(pkgconfig_path) then
            -- PKG_CONFIG_PATH typically points to .../lib/pkgconfig
            -- We want to get the prefix (two levels up) to find include and lib
            local lib_dir = path.directory(pkgconfig_path)  -- .../lib
            local prefix_dir = path.directory(lib_dir)      -- .../
            local include_dir = path.join(prefix_dir, "include")

            if os.isdir(lib_dir) then
                add_linkdirs(lib_dir)
            end
            if os.isdir(include_dir) then
                add_includedirs(include_dir)
            end
        end
    end
end

-- Options
option("examples")
    set_default(false)
    set_showmenu(true)
    set_description("Build examples")
option_end()

option("tests")
    set_default(false)
    set_showmenu(true)
    set_description("Enable tests")
option_end()

-- Define concord package (from git)
package("concord")
    add_deps("cmake")
    set_sourcedir(path.join(os.projectdir(), "build/_deps/concord-src"))

    on_fetch(function (package)
        -- Clone git repository if not exists
        local sourcedir = package:sourcedir()
        if not os.isdir(sourcedir) then
            print("Fetching concord from git...")
            os.mkdir(path.directory(sourcedir))
            os.execv("git", {"clone", "--quiet", "--depth", "1", "--branch", "2.5.0", 
                            "-c", "advice.detachedHead=false",
                            "https://github.com/onlyhead/concord.git", sourcedir})
        end
    end)

    on_install(function (package)
        local configs = {}
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" .. (package:is_debug() and "Debug" or "Release"))
        import("package.tools.cmake").install(package, configs)
    end)
package_end()

-- Add required packages
add_requires("concord")

if has_config("tests") then
    add_requires("doctest")
end

-- Main library target
target("geotiv")
    set_kind("static")

    -- Add source files
    add_files("src/*.cpp")

    -- Add header files
    add_headerfiles("include/(geotiv/**.hpp)")
    add_includedirs("include", {public = true})

    -- Link dependencies
    add_packages("concord")

    -- Set install files
    add_installfiles("include/(geotiv/**.hpp)")
    on_install(function (target)
        local installdir = target:installdir()
        os.cp(target:targetfile(), path.join(installdir, "lib", path.filename(target:targetfile())))
    end)
target_end()

-- Examples (only build when geotiv is the main project)
if has_config("examples") and os.projectdir() == os.curdir() then
    for _, filepath in ipairs(os.files("examples/*.cpp")) do
        local filename = path.basename(filepath)
        target(filename)
            set_kind("binary")
            add_files(filepath)
            add_deps("geotiv")
            add_packages("concord")
            add_includedirs("include")
        target_end()
    end
end

-- Tests (only build when geotiv is the main project)
if has_config("tests") and os.projectdir() == os.curdir() then
    for _, filepath in ipairs(os.files("test/*.cpp")) do
        local filename = path.basename(filepath)
        target(filename)
            set_kind("binary")
            add_files(filepath)
            add_deps("geotiv")
            add_packages("concord", "doctest")
            add_includedirs("include")
            add_defines("DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN")

            -- Add as test
            add_tests("default", {rundir = os.projectdir()})
        target_end()
    end
end

-- Task to generate CMakeLists.txt
task("cmake")
    on_run(function ()
        import("core.project.config")

        -- Load configuration
        config.load()

        -- Generate CMakeLists.txt
        os.exec("xmake project -k cmakelists")

        print("CMakeLists.txt generated successfully!")
    end)

    set_menu {
        usage = "xmake cmake",
        description = "Generate CMakeLists.txt from xmake.lua",
        options = {}
    }
task_end()
