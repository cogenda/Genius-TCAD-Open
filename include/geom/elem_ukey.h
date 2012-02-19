#ifndef __elem_ukey_h__
#define __elem_ukey_h__

#include <vector>
#include <algorithm>


/**
 * Unique key for a elem consisting of a variable number of vertices.
 */
class ElemKey
{
public:

  /**
   * constructor
   * build key by elem's node id
   **/
  ElemKey(const Elem * elem);

  /**
   * constructor
   * @param vertices   the vector of node id
   */
  ElemKey(const std::vector<unsigned int> &vertices);

  /**
   * constructor
   * @param elem        elem
   * @param local_ids   the vector of local node index of elem
   */
  ElemKey(const Elem * elem, const std::vector<unsigned int> &local_ids);

  /**
   * test if two ElemKey objects are equal
   */
  class Equal
  {
  public:
    /**
     * Call ElemKey::Equal()(a,b) to check equal-ness.
     * Two keys are equal iff
     *   # they have the same number of vertices
     *   # every vertex index matches
     */
    bool operator()(const ElemKey &a, const ElemKey &b) const;
  };

  /**
   * weakly less test of two ElemKey objects
   */
  class Less
  {
    public:

      /**
       * Call ElemKey::Less()(a,b) to check less.
       */
      bool operator () (const ElemKey &a, const ElemKey &b) const;
  };

  /**
   * compute the hash function of a key
   */
  class Hash
  {
  public:
    unsigned int operator()(const ElemKey &a) const;
  };

  friend std::ostream& operator<< (std::ostream& os, const ElemKey &key);

private:
  /**
   * internally stored, sorted, list of vertex indices
   */
  std::vector<unsigned int> _vertices;
};


//---------------------------------------


ElemKey::ElemKey(const Elem * elem)
{
  for(unsigned int n=0; n<elem->n_nodes(); ++n)
    _vertices.push_back(elem->node(n));
  std::sort(_vertices.begin(), _vertices.end());
}


ElemKey::ElemKey(const std::vector<unsigned int> &vertices)
{
  _vertices = vertices;
  std::sort(_vertices.begin(), _vertices.end());
}


ElemKey::ElemKey(const Elem * elem, const std::vector<unsigned int> &local_ids)
{
  for(unsigned int n=0; n<local_ids.size(); ++n)
    _vertices.push_back(elem->node(local_ids[n]));
  std::sort(_vertices.begin(), _vertices.end());
}


unsigned int ElemKey::Hash::operator()(const ElemKey &a) const
{
  /*
  One-at-a-Time hash designed by Bob Jenkins
  http://eternallyconfuzzled.com/tuts/algorithms/jsw_tut_hashing.aspx
  */

  unsigned int h = 0; // should be 32bit
  for (std::vector<unsigned int>::const_iterator it=a._vertices.begin();
       it!=a._vertices.end(); it++)
  {
    unsigned int k = *it;

    for (int i=0; i<sizeof(unsigned int); i++)
    {
      unsigned char b = k & 0xff;
      k = k >> 8;

      h += b;
      h += (h<<10);
      h ^= (h>>6);
    }
  }

  h += (h<<3);
  h ^= (h>>11);
  h += (h<<15);

  return (unsigned int)h;
}

bool ElemKey::Equal::operator()(const ElemKey &a, const ElemKey &b) const
{
  if (a._vertices.size()!=b._vertices.size()) return false;

  for (std::vector<unsigned int>::const_iterator ia=a._vertices.begin(),ib=b._vertices.begin();
       ia!=a._vertices.end() && ib!=b._vertices.end();
       ia++,ib++
      )
    if (*ia != *ib) return false;

  return true;
}

bool ElemKey::Less::operator () (const ElemKey &a, const ElemKey &b) const
{
  if (a._vertices.size() != b._vertices.size()) return (a._vertices.size() < b._vertices.size());

  for (std::vector<unsigned int>::const_iterator ia=a._vertices.begin(),ib=b._vertices.begin();
       ia!=a._vertices.end() && ib!=b._vertices.end();
       ia++,ib++
      )
    if (*ia != *ib) return (*ia) < (*ib);

  return false;
}

std::ostream& operator<< (std::ostream& os, const ElemKey &key)
{
  os << "key:" << std::hex << std::setw(10) << ElemKey::Hash()(key) << std::dec;
  os << "  vertices:";
  for (std::vector<unsigned int>::const_iterator it=key._vertices.begin();
       it!=key._vertices.end(); it++)
  {
    os << " " << *it;
  }
  return os;
}

#endif


