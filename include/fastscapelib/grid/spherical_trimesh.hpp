#ifndef FASTSCAPELIB_GRID_SPHERICAL_TRIMESH_H_
#define FASTSCAPELIB_GRID_SPHERICAL_TRIMESH_H_

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <memory>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    template <class S>
    class spherical_trimesh_xt;

    /**
     * 2-d spherical triangular mesh specialized types.
     */
    template <class S>
    struct grid_inner_types<spherical_trimesh_xt<S>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using xt_selector = S;
        static constexpr std::size_t xt_ndims = 1;

        static constexpr uint8_t n_neighbors_max = 0;
        using neighbors_cache_type = neighbors_no_cache<0>;
    };

    /**
     * @brief 2-dimensional triangular (unstructured) mesh on the Earth or a
     * planetary sphere.
     *
     * @tparam S The xtensor container selector for data array members.
     */
    template <class S>
    class spherical_trimesh_xt : public grid<spherical_trimesh_xt<S>>
    {
    public:
        using self_type = spherical_trimesh_xt<S>;
        using base_type = grid<self_type>;

        using cgal_kernel_type = typename CGAL::Exact_predicates_inexact_constructions_kernel;
        using cgal_earth_sphere_traits =
            typename CGAL::Projection_on_sphere_traits_3<cgal_kernel_type>;
        using cgal_mesh_type =
            typename CGAL::Delaunay_triangulation_on_sphere_2<cgal_earth_sphere_traits>;
        using cgal_point_3 = typename cgal_mesh_type::Point_3;

        static constexpr double earth_radius = 6.371e6;

        template <class I>
        spherical_trimesh_xt(const I& points_begin,
                             const I& points_end,
                             double radius = earth_radius);

        static spherical_trimesh_xt from_icosphere(unsigned int n_subdivisions,
                                                   double radius = earth_radius);

    protected:
        std::unique_ptr<cgal_mesh_type> m_cgal_mesh_ptr;

        friend class grid<self_type>;
    };


    /**
     * @name Constructors
     */
    //@{
    /**
     * Creates a new triangular mesh on the sphere from a given range of points.
     *
     * This basic constructor requires existing points defined as CGAL objects.
     * See other available factory methods for more convenient ways of
     * generating a new mesh.
     *
     * @param points_begin STL-compatible input iterator of CGAL ``Point_on_sphere_2``
     * or ``Point_3`` objects.
     * @param points_end STL-compatible input iterator of CGAL ``Point_on_sphere_2``
     * or ``Point_3`` objects.
     * @param radius The sphere radius (by default the approximate Earth radius in meters).
     *
     */
    template <class S>
    template <class I>
    spherical_trimesh_xt<S>::spherical_trimesh_xt(const I& points_begin,
                                                  const I& points_end,
                                                  double radius)
        : base_type(0)
    {
        cgal_earth_sphere_traits sphere_traits(cgal_point_3(0., 0., 0.), radius);
        std::make_unique<cgal_mesh_type>(points_begin, points_end, sphere_traits);
    }
    //@}

    /**
     * @name Factories
     */
    //@{
    /**
     * Creates a new triangular mesh on the sphere from an icosphere.
     *
     * The vertices of the new mesh are quasi-uniformly distributed on the sphere.
     * This works by first creating an isosahedron, then subdividing each of its
     * triangular faces into a set of smaller triangles and finally project the
     * vertices on the sphere.
     *
     * @param n_subdivisions The number of subdivisions of the initial isosahedron.
     * @param radius The sphere radius (by default the approximate Earth radius in meters).
     */
    template <class S>
    spherical_trimesh_xt<S> spherical_trimesh_xt<S>::from_icosphere(unsigned int n_subdivisions,
                                                                    double radius)
    {
        using surface_mesh_type = CGAL::Surface_mesh<cgal_point_3>;

        surface_mesh_type temp_mesh;

        CGAL::make_icosahedron<surface_mesh_type, cgal_point_3>(
            temp_mesh, cgal_point_3(0., 0., 0.), radius);
        CGAL::Subdivision_method_3::Loop_subdivision(
            temp_mesh, CGAL::parameters::number_of_iterations(n_subdivisions));

        std::vector<cgal_point_3> points;
        for (const auto& v : temp_mesh.vertices())
        {
            points.emplace_back(temp_mesh.point(v));
        }

        return spherical_trimesh_xt<S>(points.begin(), points.end(), radius);
    }
    //@}
}

#endif  // FASTSCAPELIB_GRID_SPHERICAL_TRIMESH_H_
