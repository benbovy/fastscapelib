#ifndef FASTSCAPELIB_HEALPIX_GRID_H_
#define FASTSCAPELIB_HEALPIX_GRID_H_

#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/healpix_tables.h"
#include "healpix_cxx/vec3.h"

// conflict between healpix xcomplex macro and xtl xcomplex
#undef xcomplex
#include <xtensor/xbroadcast.hpp>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/utils/xtensor_containers.hpp"

#include <math.h>
#include <memory>
#include <vector>


namespace fastscapelib
{
    namespace detail
    {
        inline double vec3_distance(vec3 a, vec3 b)
        {
            return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2)
                             + std::pow(a.z - b.z, 2));
        }

    }

    template <class S, class T>
    class healpix_grid;

    /**
     * Healpix grid specialized types
     */
    template <class S, class T>
    struct grid_inner_types<healpix_grid<S, T>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using container_selector = S;
        static constexpr std::size_t container_ndims = 1;

        static constexpr uint8_t n_neighbors_max = 8;
        using neighbors_cache_type = neighbors_no_cache<n_neighbors_max>;
    };

    /**
     * @brief 2-dimensional grid on the sphere (HEALPix).
     *
     * Fastscapelib grid adapter for a HEALPix (Hierarchical Equal Area
     * isoLatitude Pixelation of a sphere) grid.
     *
     * @tparam S The container selector for data array members.
     * @tparam T The integer type used to store the HEALPix grid node indices.
     */
    template <class S = xt_selector, class T = int>
    class healpix_grid : public grid<healpix_grid<S, T>>
    {
    public:
        using self_type = healpix_grid<S, T>;
        using base_type = grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using container_selector = typename base_type::container_selector;
        using container_type = fixed_shape_container_t<container_selector,
                                                       grid_data_type,
                                                       inner_types::container_ndims>;

        using neighbors_type = typename base_type::neighbors_type;
        using neighbors_indices_type = typename base_type::neighbors_indices_type;
        using neighbors_distances_type = typename base_type::neighbors_distances_type;

        using nodes_status_type = typename base_type::nodes_status_type;
        using nodes_status_array_type = fixed_shape_container_t<container_selector, node_status, 1>;

        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        healpix_grid(T nside,
                     const nodes_status_array_type& nodes_status,
                     double radius = numeric_constants<>::EARTH_RADIUS);

        // TODO: factory calculating nside from a given approx. cell area.

        void set_nodes_status(const nodes_status_array_type& nodes_status);

        T nside() const;
        double radius() const;

    protected:
        using healpix_type = T_Healpix_Base<T>;
        std::unique_ptr<healpix_type> m_healpix_obj_ptr;

        shape_type m_shape;
        size_type m_size;
        double m_radius;
        double m_node_area;

        nodes_status_type m_nodes_status;

        using neighbors_distances_impl_type = typename base_type::neighbors_distances_impl_type;
        using neighbors_indices_impl_type = typename base_type::neighbors_indices_impl_type;

        std::vector<size_type> m_neighbors_count;
        std::vector<neighbors_indices_impl_type> m_neighbors_indices;
        std::vector<neighbors_distances_impl_type> m_neighbors_distances;

        void set_neighbors();

        inline container_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        inline size_type neighbors_count_impl(const size_type& idx) const;

        void neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                    const size_type& idx) const;

        inline const neighbors_distances_impl_type& neighbors_distances_impl(
            const size_type& idx) const;

        static constexpr std::size_t dimension_impl() noexcept;

        friend class grid<self_type>;
    };


    /**
     * @name Constructors
     */
    /**
     * Creates a new HEALPix grid
     *
     * @param nside The number of divisions along the side of a base-resolution HEALPix pixel.
     * @param radius The radius of the sphere (default: Earth radius in meters).
     */
    template <class S, class T>
    healpix_grid<S, T>::healpix_grid(T nside,
                                     const nodes_status_array_type& nodes_status,
                                     double radius)
        : base_type(0)
        , m_radius(radius)
    {
        m_healpix_obj_ptr
            = std::make_unique<healpix_type>(nside, Healpix_Ordering_Scheme::RING, SET_NSIDE);

        m_size = static_cast<size_type>(m_healpix_obj_ptr->Npix());
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
        m_node_area = 4.0 * M_PI * m_radius * m_radius / m_size;

        // will also compute grid node neighbors
        set_nodes_status(nodes_status);
    }
    //@}

    template <class S, class T>
    void healpix_grid<S, T>::set_nodes_status(const nodes_status_array_type& nodes_status)
    {
        if (!xt::same_shape(nodes_status.shape(), m_shape))
        {
            throw std::invalid_argument(
                "invalid shape for nodes_status array (expects shape [N] where N is the total number of nodes)");
        }
        m_nodes_status = nodes_status;

        // maybe invalidates the grid node neighbors so it must be (re)computed
        set_neighbors();
    }

    template <class S, class T>
    void healpix_grid<S, T>::set_neighbors()
    {
        m_neighbors_count.resize(m_size);
        m_neighbors_indices.resize(m_size);
        m_neighbors_distances.resize(m_size);

        fix_arr<T, 8> temp_neighbors_indices;

        for (size_type inode = 0; inode < m_size; inode++)
        {
            if (m_nodes_status[inode] == node_status::ghost)
            {
                continue;
            }

            T inode_ = static_cast<T>(inode);
            size_type neighbors_count = 0;
            auto inode_vec3 = m_healpix_obj_ptr->pix2vec(inode_);

            m_healpix_obj_ptr->neighbors(inode_, temp_neighbors_indices);

            for (size_type k = 1; k < 8; ++k)
            {
                auto ineighbor = temp_neighbors_indices[k];
                auto ineighbor_ = static_cast<size_type>(ineighbor);

                if (ineighbor > -1 && m_nodes_status[ineighbor_] != node_status::ghost)
                {
                    m_neighbors_indices[inode][neighbors_count] = ineighbor_;

                    auto ineighbor_vec3 = m_healpix_obj_ptr->pix2vec(ineighbor);
                    m_neighbors_distances[inode][neighbors_count]
                        = detail::vec3_distance(inode_vec3, ineighbor_vec3);

                    neighbors_count++;
                }
            }

            m_neighbors_count[inode] = neighbors_count;
        }
    }

    template <class S, class T>
    auto healpix_grid<S, T>::nside() const -> T
    {
        return m_healpix_obj_ptr->Nside();
    }

    template <class S, class T>
    double healpix_grid<S, T>::radius() const
    {
        return m_radius;
    }

    template <class S, class T>
    inline auto healpix_grid<S, T>::nodes_areas_impl() const -> container_type
    {
        return xt::broadcast(m_node_area, m_shape);
    }

    template <class S, class T>
    inline auto healpix_grid<S, T>::nodes_areas_impl(const size_type& /*idx*/) const noexcept
        -> grid_data_type
    {
        return m_node_area;
    }

    template <class S, class T>
    inline auto healpix_grid<S, T>::neighbors_count_impl(const size_type& idx) const -> size_type
    {
        return m_neighbors_count[idx];
    }

    template <class S, class T>
    void healpix_grid<S, T>::neighbors_indices_impl(neighbors_indices_impl_type& neighbors,
                                                    const size_type& idx) const
    {
        const auto& size = m_neighbors_count[idx];

        for (size_type i = 0; i < size; i++)
        {
            neighbors[i] = m_neighbors_indices[idx][i];
        }
    }

    template <class S, class T>
    auto healpix_grid<S, T>::neighbors_distances_impl(const size_type& idx) const
        -> const neighbors_distances_impl_type&
    {
        return m_neighbors_distances[idx];
    }

    template <class S, class T>
    constexpr std::size_t healpix_grid<S, T>::dimension_impl() noexcept
    {
        return 2;
    }
}

#endif  // FASTSCAPELIB_HEALPIX_GRID_H_
