#include "gtest/gtest.h"

#include "fastscapelib/grid/spherical_trimesh.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {
        using grid_type = fs::spherical_trimesh_xt<fs::xt_selector>;

        TEST(spherical_trimesh, from_icosphere)
        {
            grid_type smesh = grid_type::from_icosphere(8);
            ASSERT_TRUE(false);
            // std::cout << smesh.cgal_mesh().number_of_vertices() << std::endl;
        }

    }

}
