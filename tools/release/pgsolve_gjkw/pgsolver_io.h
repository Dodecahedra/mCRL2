#include <limits>
#include <stdexcept>
#include <sstream>

#include "pg.h"
using namespace mcrl2::log;
namespace detail
{

template<typename ParityGame>
class pgsolver_printer
{
public:
  pgsolver_printer(ParityGame& pg) :
    m_pg(pg)
  {
  }

  void
  dump(std::ostream& s)
  {
    if (boost::num_vertices(m_pg) == 0)
      return;
    s << "parity " << num_vertices(m_pg) - 1 << ";" << std::endl;

    typename boost::graph_traits< ParityGame >::vertex_iterator i, end;
    typename boost::graph_traits< ParityGame >::adjacency_iterator ai, aend;

    for (boost::tie(i, end) = vertices(m_pg); i != end; ++i)
    {
      s << *i << ' ' << m_pg[*i].prio << ' '
        << (m_pg[*i].player == even ? '0' : '1');
      boost::tie(ai, aend) = boost::adjacent_vertices(*i, m_pg);
      if (ai != aend)
      {
        s << ' ' << *ai;
        while (++ai != aend)
        {
          s << ',' << *ai;
        }
      }
      else
        s << " \"no outgoing edges!\"";
      s << ';' << std::endl;
    }
  }

private:
  ParityGame& m_pg;
};

template<typename ParityGame>
class pgsolver_parser
{
  public:
    pgsolver_parser(ParityGame& pg) :
      m_pg(pg), m_startvertex(0)
    {
    }
    void
    load(std::istream& s)
    {
      parse_header(s);
      parse_body(s);
    }

  private:
    ParityGame& m_pg;
    size_t m_startvertex;

    void
    parse_error(std::istream& s, const char* msg)
    {
      s.clear();
      std::stringstream buf;
      std::string token;
      std::streampos p = s.tellg();
      std::streamsize l = 0;
      std::streampos b = 0;
      s >> token;
      s.seekg(0, std::ios::beg);
      while (s.tellg() < p)
      {
        ++l;
        b = s.tellg();
        s.ignore(p, '\n');
      }
      buf << msg << " Error occurred while parsing '" << token
          << "' at line " << l << ", column " << (p - b) << ".";
      throw std::runtime_error(buf.str());
    }

    void
    parse_header(std::istream& s)
    {
      std::string firstword;
      s >> std::skipws >> firstword;
      if (firstword == "parity")
      {
        size_t n;
        char c;
        s >> n;
        m_pg = ParityGame(n+1);
        s >> c;
        if (c != ';')
          parse_error(s, "Invalid header, expected semicolon.");

        // mlsolver allows start keyword
        std::string secondword;
        s >> std::skipws >> secondword;
        if (secondword == "start")
        {
          size_t start_n;
          char c;
          s >> start_n;
          if(start_n != 0)
            parse_error(s, "Invalid start vertex, expected 0.");
          s >> c;
          if (c != ';')
            parse_error(s, "Invalid header, expected semicolon.");
        }
        else
        {
          s.seekg(-secondword.length(), std::ios::cur);
          return;
        }
      }
      else
      {
        assert(firstword != "start");
        s.seekg(-firstword.length(), std::ios::cur);
        return;
      }
    }

    void
    parse_vertex(std::istream& s)
    {
      size_t index;
      size_t succ;
      char c;
      s >> index;
      if (s.fail())
      {
        if (s.eof())
          return;
        parse_error(s, "Could not parse vertex index.");
      }
      while(index >= boost::num_vertices(m_pg))
        add_vertex(m_pg);

      assert(index < boost::num_vertices(m_pg));
      typename ParityGame::vertex_descriptor v = vertex(index, m_pg);

      //m_pg[v].vertex_index = index;
      s >> m_pg[v].prio;

      if (s.fail())
        parse_error(s, "Could not parse vertex priority.");
      s >> c;
      if (s.fail() or c < '0' or c > '1')
        parse_error(s, "Could not parse vertex player.");

      m_pg[v].player = c == '0' ? even : odd;
      do
      {
        s >> succ;
        if (s.fail())
          parse_error(s, "Could not parse successor index.");
        while(succ >= boost::num_vertices(m_pg))
          boost::add_vertex(m_pg);
        typename ParityGame::vertex_descriptor u = boost::vertex(succ, m_pg);
        boost::add_edge(v, u, m_pg);

        s >> c;
        if (s.fail())
          c = ';'; // Allow missing semicolon at end of file.
      }
      while (c == ',');
      if (c == '"')
      {
        s.ignore(std::numeric_limits<std::streamsize>::max(), '"');
        s >> c;
      }
      if (c != ';')
        parse_error(s, "Invalid vertex specification, expected semicolon.");
    }

    void
    parse_body(std::istream& s)
    {
      size_t n = 0, N = boost::num_vertices(m_pg);
      N = N ? N : (size_t) -1;
      while (!s.eof() && n != N)
      {
        try
        {
          parse_vertex(s);
          ++n;
        }
        catch (std::runtime_error& e)
        {
          std::stringstream msg;
          msg << "Could not parse vertex " << n << ": " << e.what();
          throw std::runtime_error(msg.str());
        }
      }
    }
};

} // namespace detail

template<typename ParityGame>
inline void
parse_pgsolver(ParityGame& graph, std::istream& s, mcrl2::utilities::execution_timer& timer)
{
  mCRL2log(verbose)
    << "Loading parity game." << std::endl;
  timer.start("load");
  detail::pgsolver_parser<ParityGame> parser(graph);
  parser.load(s);
  timer.finish("load");
  mCRL2log(verbose)
    << "Parity game contains " << boost::num_vertices(graph) << " nodes and "
        << boost::num_edges(graph) << " edges." << std::endl;
}

template<typename ParityGame>
inline void
print_pgsolver(ParityGame& graph, std::ostream& os)
{
  mCRL2log(verbose) << "Printing parity game." << std::endl;
  detail::pgsolver_printer<ParityGame> p(graph);
  p.dump(os);
}