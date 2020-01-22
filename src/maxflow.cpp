//
// Created by Jan Groschaft on 8.3.19.
//


#include <ios>
#include <fstream>
#include "common_types.h"
#include "command_line_parser.h"
#include "graph_loader.h"
#include "measure.h"
#include "algorithms/parallel/push_relabel_segment.h"
#include "algorithms/sequential/ahuja_orlin.h"
#include "algorithms/parallel/parallel_push_relabel.h"
#include "algorithms/sequential/push_relabel_highest.h"
#include "algorithms/sequential/push_relabel_fifo.h"
#include "algorithms/sequential/edmonds_karp.h"
#include "algorithms/sequential/dinic.h"
#include "algorithms/parallel/ahuja_orlin_segment.h"


std::istream & get_input_stream ( std::ifstream & in_file, std::optional<std::string> & file_path )
{
    if ( file_path . has_value () )
    {
        in_file . open ( *file_path );
        return in_file;
    }
    return std::cin;
}


template <typename U>
void print_result ( const measurement_result<U> & res, std::string_view solver, std::string_view filename,
        std::size_t thr_cnt )
{
    std::cout << "solver:\t\t" << solver << '\n' <<
              "filename:\t" << filename << '\n' <<
              "flow:\t\t" << res . max_flow << '\n' <<
              "time read:\t" << res . time_read . count () << " ms\n" <<
              "time init:\t" << res . time_init . count () << " ms\n" <<
              "time solve:\t" << res . time_solve . count () << " ms\n" <<
              "# of threads:\t" << thr_cnt << "\n";
}

auto load_graph_and_run ( std::istream & is, solver solver, std::size_t & thread_count )
{
    using T = uint32_t;
    using U = uint64_t;

    auto get_graph = [] ( std::istream & is )
    {
        auto start = std::chrono::high_resolution_clock::now ();
        auto[graph, source, sink] = load_graph<T, U, basic_edge> ( is );
        auto end = std::chrono::high_resolution_clock::now ();
        return std::make_tuple ( graph, source, sink, std::chrono::duration_cast<std::chrono::milliseconds> ( end - start ) );
    };

    auto get_graph_with_cached_edge = [] ( std::istream & is )
    {
        auto start = std::chrono::high_resolution_clock::now ();
        auto[graph, source, sink] = load_graph<T, U, cached_edge> ( is );
        set_reverse_edge_cap ( graph );
        auto end = std::chrono::high_resolution_clock::now ();
        return std::make_tuple ( graph, source, sink, std::chrono::duration_cast<std::chrono::milliseconds> ( end - start ) );
    };

    measurement_result<U> result {};

    switch ( solver )
    {
        case solver::ek:
        {
            auto[graph, source, sink, time_read] = get_graph ( is );
            result = measure_sequential<edmonds_karp::max_flow_instance> ( std::move ( graph ), source, sink );
            result . time_read = time_read;
            thread_count = 1;
            break;
        }
        case solver::din:
        {
            auto[graph, source, sink, time_read] = get_graph ( is );
            result = measure_sequential<dinic::max_flow_instance> ( std::move ( graph ), source, sink );
            result . time_read = time_read;
            thread_count = 1;
            break;
        }

        case solver::prf:
        {
            auto[graph, source, sink, time_read] = get_graph_with_cached_edge ( is );
            result = measure_sequential<push_relabel_fifo::max_flow_instance> ( std::move ( graph ), source, sink );
            result . time_read = time_read;
            thread_count = 1;
            break;
        }
        case solver::prh:
        {
            auto[graph, source, sink, time_read] = get_graph_with_cached_edge ( is );
            result = measure_sequential<push_relabel_highest::max_flow_instance> ( std::move ( graph ), source, sink );
            result . time_read = time_read;
            thread_count = 1;
            break;
        }
        case solver::ppr:
        {
            auto[graph, source, sink, time_read] = get_graph_with_cached_edge ( is );
            result = measure_parallel<parallel_push_relabel::max_flow_instance> ( std::move ( graph ), source, sink, thread_count );
            result . time_read = time_read;
            break;
        }
        case solver::prs:
        {
            auto[graph, source, sink, time_read] = get_graph_with_cached_edge ( is );
            result = measure_parallel<push_relabel_segment::max_flow_instance> ( std::move ( graph ), source, sink, thread_count );
            result . time_read = time_read;
            break;
        }
        case solver::ao:
        {
            auto[graph, source, sink, time_read] = get_graph_with_cached_edge ( is );
            result = measure_sequential<ahuja_orlin::max_flow_instance> ( std::move ( graph ), source, sink );
            result . time_read = time_read;
            thread_count = 1;
            break;
        }
        case solver::aos:
        {
            auto[graph, source, sink, time_read] = get_graph_with_cached_edge ( is );
            result = measure_parallel<ahuja_orlin_segment::max_flow_instance> ( std::move ( graph ), source, sink, thread_count );
            result . time_read = time_read;
            break;
        }
        default:
            throw std::logic_error ( "Unknown solver" );
    }
    return result;
}


int main ( int argc, char * argv[] )
{
    std::ios_base::sync_with_stdio ( false );
    command_line_parser parser;
    if ( !parser . parse_arguments ( argc, argv ) )
        return 1;
    auto solver = parser . get_solver ();
    auto solver_str = parser . get_solver_str ();
    auto file_path = parser . get_filename ();
    auto thr_cnt = parser . get_thread_count ();

    std::ifstream in_file;
    auto & stream = get_input_stream ( in_file, file_path );
    if ( !stream )
    {
        std::cerr << "Unable to open file " << *file_path << '\n';
        return 1;
    }

    //get filename without full path, or stdin
    std::string filename = "stdin";
    if ( file_path . has_value () )
    {
        filename = (*file_path) . substr((*file_path).find_last_of("/\\") + 1);
    }

    //run
    auto result = load_graph_and_run ( stream, solver, thr_cnt );
    print_result ( result, solver_str, filename, thr_cnt );
    return 0;
}
