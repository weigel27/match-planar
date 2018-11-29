// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/lattice/bgl_properties.h,v 1.3 2008-08-11 12:19:46 weigel Exp $ 

#ifndef _BGL_PROPERTIES_H_
#define _BGL_PROPERTIES_H_

#include <boost/graph/properties.hpp>

enum vertex_content_t { vertex_content };
namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, content);
}

enum edge_content_t { edge_content };
namespace boost {
  BOOST_INSTALL_PROPERTY(edge, content);
}

enum edge_dual_t { edge_dual };
namespace boost {
  BOOST_INSTALL_PROPERTY(edge, dual);
}

enum edge_predecessortree_t { edge_predecessortree };
namespace boost {
  BOOST_INSTALL_PROPERTY(edge, predecessortree);
}

#endif
