#include <fstream>
#include <stack>

#include "octree.h"



OcTree::OcTree(const Point &b1, const Point &b2, OcTreeDataBase * data, int depth)
    :tree<OcTreeNode >()
{
  const Point * p0 = this->add_point ( Point(b1.x(), b1.y(), b1.z()) );
  const Point * p1 = this->add_point ( Point(b2.x(), b1.y(), b1.z()) );
  const Point * p2 = this->add_point ( Point(b2.x(), b2.y(), b1.z()) );
  const Point * p3 = this->add_point ( Point(b1.x(), b2.y(), b1.z()) );
  const Point * p4 = this->add_point ( Point(b1.x(), b1.y(), b2.z()) );
  const Point * p5 = this->add_point ( Point(b2.x(), b1.y(), b2.z()) );
  const Point * p6 = this->add_point ( Point(b2.x(), b2.y(), b2.z()) );
  const Point * p7 = this->add_point ( Point(b1.x(), b2.y(), b2.z()) );

  OcTreeNode root_node;
  root_node.set_corner_point(p0, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  root_node.set_corner_point(p1, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  root_node.set_corner_point(p2, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  root_node.set_corner_point(p3, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  root_node.set_corner_point(p4, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  root_node.set_corner_point(p5, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  root_node.set_corner_point(p6, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  root_node.set_corner_point(p7, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  root_node.set_data(data);

  _data.insert(data);

  this->set_head ( root_node );

  for(int d=0; d<depth; ++d)
  {
    tree_leaf_iterator leaf_it = tree< OcTreeNode >::begin_leaf();
    for ( ; leaf_it != tree< OcTreeNode >::end_leaf(); )
    {
      tree_iterator_base this_leaf = leaf_it++;
      subdivide ( this_leaf );
    }
  }
}


OcTree::OcTree ( const OcTreeNode & root_node )
    :tree<OcTreeNode > ( root_node )
{
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Bottom, OcTreeLocation::L_Back)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Top,    OcTreeLocation::L_Back)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Bottom, OcTreeLocation::L_Front)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Top,    OcTreeLocation::L_Front)) );
  this->insert_point( root_node.get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front)) );
}



OcTree::~OcTree()
{
  for ( unsigned int n=0; n<_points.size(); ++n )
    delete _points[n];

  std::set<OcTreeDataBase *>::const_iterator it = _data.begin();
  for(; it != _data.end(); ++it)
    delete *it;
}



const Point * OcTree::add_point ( const Point & point )
{
  //for(unsigned int n=0; n<_points.size(); ++n)
  //  if( point.absolute_fuzzy_equals(*_points[n]) )
  //    return _points[n];
  if ( _point_to_id.find ( &point ) != _point_to_id.end() )
    return _point_to_id.find ( &point )->first;

  Point * new_point = new Point ( point );
  _point_to_id[new_point] = _points.size();
  _points.push_back ( new_point );
  return new_point;
}


void OcTree::insert_point( const Point * p )
{
  _point_to_id[p] = _points.size();
  _points.push_back ( p );
}



void OcTree::subdivide ( tree_iterator_base & leaf_it )
{
  /*
    Leaf: 7        6
          o--------o
         /:       /|
        / :      / |
     4 /  :   5 /  |
      o--------o   |
      |   o....|...o 2
      |  .3    |  /
      | .      | /
      |.       |/
      o--------o
      0        1
  */
  const Point * p0 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  const Point * p1 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  const Point * p2 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  const Point * p3 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  const Point * p4 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  const Point * p5 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  const Point * p6 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  const Point * p7 = leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));


  const Point *c  = this->add_point ( ( *p0 + *p6 ) *0.5 );
  const Point *e0 = this->add_point ((*p0 + *p1)*0.5);
  const Point *e1 = this->add_point ((*p1 + *p2)*0.5);
  const Point *e2 = this->add_point ((*p2 + *p3)*0.5);
  const Point *e3 = this->add_point ((*p0 + *p3)*0.5);

  const Point *e4 = this->add_point ((*p0 + *p4)*0.5);
  const Point *e5 = this->add_point ((*p1 + *p5)*0.5);
  const Point *e6 = this->add_point ((*p2 + *p6)*0.5);
  const Point *e7 = this->add_point ((*p3 + *p7)*0.5);

  const Point *e8 = this->add_point ((*p4 + *p5)*0.5);
  const Point *e9 = this->add_point ((*p5 + *p6)*0.5);
  const Point *e10 = this->add_point ((*p6 + *p7)*0.5);
  const Point *e11 = this->add_point ((*p7 + *p4)*0.5);

  const Point *f0 = this->add_point ((*p0 + *p1 + *p2 + *p3)*0.25);
  const Point *f1 = this->add_point ((*p0 + *p1 + *p5 + *p4)*0.25);
  const Point *f2 = this->add_point ((*p1 + *p2 + *p6 + *p5)*0.25);
  const Point *f3 = this->add_point ((*p3 + *p2 + *p6 + *p7)*0.25);
  const Point *f4 = this->add_point ((*p0 + *p3 + *p7 + *p4)*0.25);
  const Point *f5 = this->add_point ((*p4 + *p5 + *p6 + *p7)*0.25);


  OcTreeNode child0(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child0.set_corner_point(p0, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child0.set_corner_point(e0, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child0.set_corner_point(f0, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child0.set_corner_point(e3, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child0.set_corner_point(e4, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child0.set_corner_point(f1, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child0.set_corner_point(c,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child0.set_corner_point(f4, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child0 = append_child ( leaf_it,  child0);

  OcTreeNode child1(OcTreeLocation(OcTreeLocation::L_Right,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child1.set_corner_point(e0, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child1.set_corner_point(p1, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child1.set_corner_point(e1, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child1.set_corner_point(f0, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child1.set_corner_point(f1, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child1.set_corner_point(e5, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child1.set_corner_point(f2, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child1.set_corner_point(c,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child1 = append_child ( leaf_it,  child1);

  OcTreeNode child2(OcTreeLocation(OcTreeLocation::L_Right,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child2.set_corner_point(f0, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child2.set_corner_point(e1, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child2.set_corner_point(p2, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child2.set_corner_point(e2, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child2.set_corner_point(c,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child2.set_corner_point(f2, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child2.set_corner_point(e6, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child2.set_corner_point(f3, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child2 = append_child ( leaf_it,  child2);

  OcTreeNode child3(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child3.set_corner_point(e3, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child3.set_corner_point(f0, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child3.set_corner_point(e2, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child3.set_corner_point(p3, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child3.set_corner_point(f4, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child3.set_corner_point(c,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child3.set_corner_point(f3, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child3.set_corner_point(e7, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child3 = append_child ( leaf_it,  child3);


  OcTreeNode child4(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top, OcTreeLocation::L_Front));
  child4.set_corner_point(e4,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child4.set_corner_point(f1,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child4.set_corner_point(c,   OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child4.set_corner_point(f4,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child4.set_corner_point(p4,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child4.set_corner_point(e8,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child4.set_corner_point(f5,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child4.set_corner_point(e11, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child4 = append_child ( leaf_it,  child4);

  OcTreeNode child5(OcTreeLocation(OcTreeLocation::L_Right,   OcTreeLocation::L_Top, OcTreeLocation::L_Front));
  child5.set_corner_point(f1,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child5.set_corner_point(e5,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child5.set_corner_point(f2,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child5.set_corner_point(c,   OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child5.set_corner_point(e8,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child5.set_corner_point(p5,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child5.set_corner_point(e9,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child5.set_corner_point(f5,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child5 = append_child ( leaf_it,  child5);

  OcTreeNode child6(OcTreeLocation(OcTreeLocation::L_Right,   OcTreeLocation::L_Top, OcTreeLocation::L_Back));
  child6.set_corner_point(c,   OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child6.set_corner_point(f2,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child6.set_corner_point(e6,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child6.set_corner_point(f3,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child6.set_corner_point(f5,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child6.set_corner_point(e9,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child6.set_corner_point(p6,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child6.set_corner_point(e10, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child6 = append_child ( leaf_it,  child6);

  OcTreeNode child7(OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top, OcTreeLocation::L_Back));
  child7.set_corner_point(f4,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child7.set_corner_point(c,   OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front));
  child7.set_corner_point(f3,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child7.set_corner_point(e7,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Bottom, OcTreeLocation::L_Back));
  child7.set_corner_point(e11, OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child7.set_corner_point(f5,  OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front));
  child7.set_corner_point(e10, OcTreeLocation(OcTreeLocation::L_Right,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  child7.set_corner_point(p7,  OcTreeLocation(OcTreeLocation::L_Left,   OcTreeLocation::L_Top,    OcTreeLocation::L_Back));
  tree_iterator_base  it_child7 = append_child ( leaf_it,  child7);

  if(leaf_it->data())
  {
    std::vector<OcTreeDataBase *> new_leaf_data = leaf_it->data()->subdivide(*leaf_it);
    it_child0->set_data(new_leaf_data[0]);
    it_child1->set_data(new_leaf_data[1]);
    it_child2->set_data(new_leaf_data[2]);
    it_child3->set_data(new_leaf_data[3]);
    it_child4->set_data(new_leaf_data[4]);
    it_child5->set_data(new_leaf_data[5]);
    it_child6->set_data(new_leaf_data[6]);
    it_child7->set_data(new_leaf_data[7]);

    for(unsigned int n=0; n<new_leaf_data.size(); ++n)
      _data.insert(new_leaf_data[n]);

    _data.erase(leaf_it->data());
    delete leaf_it->data();
    leaf_it->set_data(0);
  }

}




void OcTree::refine()
{
  unsigned int flag;

  do
  {
    flag = 0;
    tree_leaf_iterator leaf_it = tree< OcTreeNode >::begin_leaf();
    for ( ; leaf_it != tree< OcTreeNode >::end_leaf(); ++leaf_it )
    {
      tree_iterator_base leaf = leaf_it;

      std::vector<OcTreeNode> neighbors;
      {
        tree_iterator_base n_top = find_neighbor(leaf_it, OcTreeLocation::L_Top);
        if(n_top.node) neighbors.push_back(*n_top);

        tree_iterator_base n_bot = find_neighbor(leaf_it, OcTreeLocation::L_Bottom);
        if(n_bot.node) neighbors.push_back(*n_bot);

        tree_iterator_base n_left = find_neighbor(leaf_it, OcTreeLocation::L_Left);
        if(n_left.node) neighbors.push_back(*n_left);

        tree_iterator_base n_right = find_neighbor(leaf_it, OcTreeLocation::L_Right);
        if(n_right.node) neighbors.push_back(*n_right);

        tree_iterator_base n_front = find_neighbor(leaf_it, OcTreeLocation::L_Front);
        if(n_front.node) neighbors.push_back(*n_front);

        tree_iterator_base n_back = find_neighbor(leaf_it, OcTreeLocation::L_Back);
        if(n_back.node) neighbors.push_back(*n_back);
      }

      OcTreeDataBase * data = leaf->data();
      if(data && data->refine(*leaf, neighbors))
      {
        leaf->divide_flag() = true;
      }
    }

    for ( leaf_it = tree< OcTreeNode >::begin_leaf(); leaf_it != tree< OcTreeNode >::end_leaf(); )
    {
      tree_iterator_base this_leaf = leaf_it++;
      if ( this_leaf->divide_flag() )
      {
        subdivide ( this_leaf );
        flag++;
        this_leaf->divide_flag() = false;
      }
    }
  }
  while ( flag );

  this->balance();

}




void OcTree::balance()
{
  unsigned int flag;

  do
  {
    flag = 0;
    tree_leaf_iterator leaf_it = tree< OcTreeNode >::begin_leaf();
    for ( ; leaf_it != tree< OcTreeNode >::end_leaf(); ++leaf_it )
    {
      tree_iterator_base leaf = leaf_it;
      int leaf_depth = depth ( leaf );

      std::vector<tree_iterator_base> neighbors;

      if ( leaf_it->get_location().has_location ( OcTreeLocation::L_Top ) )
        neighbors.push_back ( find_neighbor ( leaf, OcTreeLocation::L_Top ) );
      if ( leaf_it->get_location().has_location ( OcTreeLocation::L_Bottom ) )
        neighbors.push_back ( find_neighbor ( leaf, OcTreeLocation::L_Bottom ) );
      if ( leaf_it->get_location().has_location ( OcTreeLocation::L_Left ) )
        neighbors.push_back ( find_neighbor ( leaf, OcTreeLocation::L_Left ) );
      if ( leaf_it->get_location().has_location ( OcTreeLocation::L_Right ) )
        neighbors.push_back ( find_neighbor ( leaf, OcTreeLocation::L_Right ) );
      if ( leaf_it->get_location().has_location ( OcTreeLocation::L_Front ) )
        neighbors.push_back ( find_neighbor ( leaf, OcTreeLocation::L_Front ) );
      if ( leaf_it->get_location().has_location ( OcTreeLocation::L_Back ) )
        neighbors.push_back ( find_neighbor ( leaf, OcTreeLocation::L_Back ) );

      for ( unsigned int n=0; n<neighbors.size(); ++n )
      {
        if ( !neighbors[n].node ) continue;
        int neighbor_depth = depth ( neighbors[n] );
        if ( abs ( leaf_depth-neighbor_depth ) >1 )
        {
          neighbors[n]->divide_flag() = true;
        }
      }
    }

    for ( leaf_it = tree< OcTreeNode >::begin_leaf(); leaf_it != tree< OcTreeNode >::end_leaf(); )
    {
      tree_iterator_base this_leaf = leaf_it++;
      if ( this_leaf->divide_flag() )
      {
        subdivide ( this_leaf );
        flag++;
        this_leaf->divide_flag() = false;
      }
    }

  }
  while ( flag );

}






OcTree::tree_iterator_base OcTree::find_neighbor ( const tree_iterator_base & it, OcTreeLocation::Location x ) const
{
  std::stack<OcTreeLocation> path;
  tree_iterator_base q = it;

  while ( !is_root ( q ) )
  {
    const OcTreeLocation & loc = q->get_location();
    path.push ( loc );
    q=parent ( q );
    if ( !loc.has_location ( x ) ) break;
    if ( is_root ( q ) && loc.has_location ( x ) ) return tree_iterator_base ( 0 );
  };

  while ( !path.empty() )
  {
    OcTreeLocation loc = path.top();
    path.pop();
    if ( x==OcTreeLocation::L_Left || x==OcTreeLocation::L_Right )
      q = goto_octree_child ( q, OcTreeLocation::x_symmetry ( loc ) );
    if ( x==OcTreeLocation::L_Top || x==OcTreeLocation::L_Bottom )
      q = goto_octree_child ( q, OcTreeLocation::y_symmetry ( loc ) );
    if ( x==OcTreeLocation::L_Front || x==OcTreeLocation::L_Back )
      q = goto_octree_child ( q, OcTreeLocation::z_symmetry ( loc ) );
    if ( q.node==0 || is_child ( q ) ) return q;
  }

  return q;
}


OcTree::tree_iterator_base OcTree::goto_octree_child ( const tree_iterator_base & it,
    const OcTreeLocation & location ) const
{
  if ( it.node==0 ) return  tree_iterator_base ( 0 );
  for ( tree_sibling_iterator child_it= begin ( it ); child_it!=end ( it ); ++child_it )
    if ( child_it->get_location() ==location )
      return tree_iterator_base ( child_it );
  return  tree_iterator_base ( 0 );
}



OcTree::tree_iterator_base OcTree::find_leaf_has_point ( const Point & p ) const
{
  tree_iterator_base q(tree<OcTreeNode>::head->next_sibling);

  //std::cout<<*q;
  //std::cout<<p;
  if(!q->has_point(p)) return tree_iterator_base(0);

  while( !is_child(q) )
  {
    tree_sibling_iterator sit = tree<OcTreeNode>::begin(q);
    tree_sibling_iterator sit_end = tree<OcTreeNode>::end(q);
    for( ; sit != sit_end; ++sit )
      if( sit->has_point(p) )
    {
      q = sit;
      break;
    }
  }
  return q;
}



void OcTree::intersect(const Point &p1, const Point &p2, std::vector<std::pair<OcTreeNode, double> > & result)
{
  Point begin = p1;
  Point end   = p2;
  const Point dir = (p2-p1).unit();
  const double length_orin = (p2-p1).size();
  const double eps = 1e-6*length_orin;


  // cut edge (p1, p2) by root
  tree_iterator_base root(tree<OcTreeNode>::head->next_sibling);
  //begin point not inside root
  if(!root->has_point(begin))
  {
    std::pair<double, double> t;
    if(!root->hit(begin, dir, t)) return; // no hit point
    if(t.first > length_orin) return; // no hit point
    begin = begin + dir*t.first;

    
    //std::cout<<t.first << " " << t.second << std::endl;
    //std::cout<<"dir " << dir;
    //std::cout<<"begin " << begin;

    if(!root->has_point(begin+ dir*eps)) return; // safe guard
  }
  //end point not inside root
  if(!root->has_point(end))
  {
    std::pair<double, double> t;
    root->hit(begin+dir*eps, dir, t);
    end = begin + dir*(t.second-eps);
  }


  double length = (end-begin).size();
  // start with begin point
  do
  {
    begin = begin + dir*eps;
    length = (end-begin).size();

    tree_iterator_base q = find_leaf_has_point(begin);
    if(!q.node) return; // error exit

    std::pair<double, double> t;
    if(!q->hit(begin, dir, t)) return; //error exit

    double l = std::min(length, t.second);
    result.push_back( std::make_pair(*q, l) );

    length -= l;
    begin = begin + dir*l;

  }while(length > 0.0);

}



std::string OcTree::key( const tree_iterator_base & leaf_it) const
{
  std::stack<OcTreeLocation> path;
  tree_iterator_base q = leaf_it;

  if ( q.node==0 ) return "root";

  while ( !is_root ( q ) )
  {
    const OcTreeLocation &loc = q->get_location();
    path.push ( loc );
    q=parent ( q );
  };

  std::stringstream ss;

  while ( !path.empty() )
  {
    OcTreeLocation loc = path.top();
    path.pop();
    ss << loc.key();
  }
  return ss.str();
}



void OcTree::print_path ( const tree_iterator_base & it ) const
{
  std::stack<OcTreeLocation> path;
  tree_iterator_base q = it;

  if ( q.node==0 ) return;

  while ( !is_root ( q ) )
  {
    const OcTreeLocation &loc = q->get_location();
    path.push ( loc );
    q=parent ( q );
  };

  std::cout<<'/';

  while ( !path.empty() )
  {
    OcTreeLocation loc = path.top();
    path.pop();
    loc.print ( std::cout );
    std::cout<<'/';
  }

  std::cout<<std::endl<<std::endl;
}



void OcTree::export_vtk ( const std::string & file )
{

  std::ofstream fout;
  fout.open ( file.c_str(), std::ofstream::trunc );

  fout << "# vtk DataFile Version 3.0" <<'\n';
  fout << "Date calculated by OCTREE"    <<'\n';
  fout << "ASCII"                      <<'\n';
  fout << "DATASET UNSTRUCTURED_GRID"  <<'\n';
  fout << "POINTS " << _points.size()  << " float" << '\n';

  for ( unsigned int i=0; i<_points.size(); ++i )
  {
    fout << _points[i]->x() << " \t" << _points[i]->y()  << " \t " << _points[i]->z() << '\n';
  }

  fout << std::endl;

  unsigned int n_leafs=0;
  tree_leaf_iterator leaf_it = tree< OcTreeNode >::begin_leaf();
  for ( ; leaf_it != tree< OcTreeNode >::end_leaf(); ++leaf_it )
    n_leafs++;

  fout<<"CELLS "<<n_leafs<<" "<<9*n_leafs<<'\n';

  for ( leaf_it = tree< OcTreeNode >::begin_leaf(); leaf_it != tree< OcTreeNode >::end_leaf(); ++leaf_it )
  {
    fout << 8 << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Back)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Bottom, OcTreeLocation::L_Back)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Top,    OcTreeLocation::L_Back)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Top,    OcTreeLocation::L_Back)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Bottom, OcTreeLocation::L_Front)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Bottom, OcTreeLocation::L_Front)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Right, OcTreeLocation::L_Top,    OcTreeLocation::L_Front)) ] << " "
    << _point_to_id[leaf_it->get_point(OcTreeLocation(OcTreeLocation::L_Left,  OcTreeLocation::L_Top,    OcTreeLocation::L_Front)) ] << "\n";
  }

  fout << std::endl;

  fout << "CELL_TYPES " << n_leafs << '\n';

  for ( unsigned int i=0; i<n_leafs; ++i )
    fout << 12 << std::endl; //VTK_HEXAHEDRON
  fout << std::endl;

  fout<<"CELL_DATA " << n_leafs      <<'\n';
  fout << "SCALARS "<<"data"<<" float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;

  for ( leaf_it = tree< OcTreeNode >::begin_leaf(); leaf_it != tree< OcTreeNode >::end_leaf(); ++leaf_it )
  {
    if(leaf_it->data())
      fout <<leaf_it->data()->value("data")/leaf_it->volume() << "\n";
    else
      fout << 0.0 << "\n";
  }

  fout.close();


}

