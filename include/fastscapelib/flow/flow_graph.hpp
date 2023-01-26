/**
 * @brief Class used to compute or follow flow routes on
 * the topographic surface.
 *
 * It provides a single API to chain flow_router and
 * sink_resolver methods and access their results.
 *
 */
#ifndef FASTSCAPELIB_FLOW_FLOW_GRAPH_H
#define FASTSCAPELIB_FLOW_FLOW_GRAPH_H

#include <map>
#include <stdexcept>
#include <type_traits>
#include <string>

#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    /**
     * Main class used to compute or follow flow routes on
     * the topographic surface.
     *
     * @tparam G The grid type
     * @tparam Tag The flow graph implementation tag
     * @tparam S The xtensor selector type
     */
    template <class G, class S = typename G::xt_selector, class Tag = flow_graph_fixed_array_tag>
    class flow_graph
    {
    public:
        using self_type = flow_graph<G, S, Tag>;
        using grid_type = G;
        using xt_selector = S;
        using impl_type = detail::flow_graph_impl<G, S, Tag>;
        using operators_type = flow_operator_sequence<impl_type>;

        using size_type = typename grid_type::size_type;

        using data_type = typename grid_type::grid_data_type;
        using data_array_type = xt_array_t<xt_selector, data_type>;
        using shape_type = typename data_array_type::shape_type;
        using data_array_size_type = xt_array_t<xt_selector, size_type>;

        using graph_map = std::map<std::string, self_type>;
        using graph_impl_map = std::map<std::string, impl_type&>;
        using elevation_map = std::map<std::string, data_array_type>;

        flow_graph(G& grid, operators_type operators)
            : m_grid(grid)
            , m_impl(grid)
            , m_operators(std::move(operators))
        {
            // sanity checks
            if (!m_operators.graph_updated())
            {
                throw std::invalid_argument(
                    "must have at least one operator that updates the flow graph");
            }
            if (m_operators.out_flowdir() == flow_direction::undefined)
            {
                throw std::invalid_argument(
                    "must have at least one operator that defines the output flow direction type");
            }

            // pre-allocate graph and elevation snapshots
            for (const auto& key : operators.graph_snapshot_keys())
            {
                // m_graph_snapshots.insert({ key, self_type(grid) });
                // m_graph_impl_snapshots.insert({ key, m_graph_snapshots.at(key).m_impl });
            }
            for (const auto& key : operators.elevation_snapshot_keys())
            {
                // m_elevation_snapshots.insert({ key, {} });
            }

            // pre-allocate hydrologically corrected elevation
            if (m_operators.elevation_updated())
            {
                // TODO:
            }
        }

        const data_array_type& update_routes(const data_array_type& elevation)
        {
            data_array_type& elevation_copy = elevation;

            // reset hydrologically corrected elevation
            if (m_operators.elevation_updated())
            {
                // TODO: copy elevation values
            }

            for (const auto& op : m_operators)
            {
                op.apply(m_impl, elevation_copy);
                op.save(m_impl, m_graph_impl_snapshots, elevation_copy, m_elevation_snapshots);
            }

            return elevation;
        }

        grid_type& grid() const
        {
            return m_grid;
        }

        size_type size() const
        {
            return m_grid.size();
        }

        shape_type grid_shape() const
        {
            // grid shape may have a different type (e.g., from xtensor containers)
            auto shape = m_grid.shape();
            shape_type data_array_shape(shape.begin(), shape.end());
            return data_array_shape;
        }

        const impl_type& impl() const
        {
            return m_impl;
        }

        void accumulate(data_array_type& acc, const data_array_type& src) const
        {
            return m_impl.accumulate(acc, src);
        }
        void accumulate(data_array_type& acc, data_type src) const
        {
            return m_impl.accumulate(acc, src);
        }
        data_array_type accumulate(const data_array_type& src) const
        {
            return m_impl.accumulate(src);
        }
        data_array_type accumulate(data_type src) const
        {
            return m_impl.accumulate(src);
        }

        data_array_size_type basins()
        {
            data_array_size_type basins = data_array_type::from_shape(m_grid.shape());
            auto basins_flat = xt::flatten(basins);

            m_impl.compute_basins();
            basins_flat = m_impl.basins();

            return basins;
        }

        // used internally for creating graph snapshots
        // TODO: make it private
        flow_graph(grid_type& grid)
            : m_grid(grid)
            , m_impl(grid)
        {
        }

    private:
        grid_type& m_grid;
        impl_type m_impl;
        data_array_type m_hydro_elevation;

        graph_map m_graph_snapshots;
        graph_impl_map m_graph_impl_snapshots;
        elevation_map m_elevation_snapshots;

        operators_type m_operators;
    };


    template <class G,
              class O,
              class S = typename G::xt_selector,
              class Tag = flow_graph_fixed_array_tag>
    flow_graph<G, S, Tag> make_flow_graph(G& grid, O operators)
    {
        return flow_graph<G, S, Tag>(grid, operators);
    }
}

#endif
