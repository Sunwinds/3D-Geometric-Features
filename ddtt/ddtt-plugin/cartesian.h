/***************************************************************************
 * blitz/array/cartesian.h  Cartesian product of indirection containers
 *
 * Copyright (C) 1997-2001 Todd Veldhuizen <tveldhui@oonumerics.org>
 *
 * This code was relicensed under the modified BSD license for use in SciPy
 * by Todd Veldhuizen (see LICENSE.txt in the weave directory).
 *
 *
 * Suggestions:          blitz-dev@oonumerics.org
 * Bugs:                 blitz-bugs@oonumerics.org
 *
 * For more information, please see the Blitz++ Home Page:
 *    http://oonumerics.org/blitz/
 *
 ****************************************************************************/
#ifndef BZ_ARRAY_CARTESIAN_H
#define BZ_ARRAY_CARTESIAN_H

namespace blitz{

#define _bz_typename typename

/*
 * CartesianProduct<T_tuple,T_container> is an adaptor which represents
 * the cartesian product of several containers.  
 */

// forward declaration of iterator
template<typename T_tuple, typename T_container>
class CartesianProductIterator;

struct _cp_end_tag { };

template<typename T_tuple, typename T_container>
class CartesianProduct {
public:
    typedef T_tuple value_type;
    typedef T_tuple& reference;
    typedef const T_tuple& const_reference;
    typedef CartesianProductIterator<T_tuple,T_container> iterator;
    typedef int difference_type;
    typedef int size_type;

    iterator begin()
    { return iterator(*this, N_containers); }

    iterator end()
    { return iterator(_cp_end_tag()); }

	template<typename VECTOR>
	CartesianProduct(const VECTOR& containers)
	{
		N_containers = containers.size();

		containers_.resize(N_containers, NULL);

		for(size_t i = 0; i < containers.size(); i++)
			containers_[i] = &containers[i];
	}

	CartesianProduct(){}

    const T_container& operator[](int i)
    { return *(containers_[i]); }

protected:
    std::vector<const T_container*> containers_; 
	size_t	N_containers;
};

template<typename T_tuple, typename T_container>
class CartesianProductIterator {
public:
    typedef _bz_typename T_container::const_iterator citerator;
    typedef CartesianProductIterator<T_tuple,T_container> iterator;
    typedef CartesianProduct<T_tuple,T_container> T_cp;

    CartesianProductIterator(T_cp& container, size_t N)
    {
		N_containers = N;

		iters_.resize(N);
		firstiters_.resize(N);
		enditers_.resize(N);

        for (int i=0; i < N_containers; ++i)
        {
            firstiters_[i] = container[i].begin();
            iters_[i] = firstiters_[i];
            enditers_[i] = container[i].end();

            tuple_.push_back( *iters_[i] );
        }

        endflag_ = false;
    }

    void operator++();

    CartesianProductIterator(_cp_end_tag)
    {
        endflag_ = true;
    }

    bool operator==(const iterator& x) const
    {
        return (endflag_ == x.endflag_);
    }

    bool operator!=(const iterator& x) const
    {   
        return endflag_ != x.endflag_;
    }

    const T_tuple& operator*() const
    { return tuple_; }

protected:
    std::vector<citerator> iters_;
    std::vector<citerator> firstiters_;
    std::vector<citerator> enditers_;
    T_tuple   tuple_;
    bool      endflag_;
	int		  N_containers;
};

template<typename T_tuple, typename T_container>
void CartesianProductIterator<T_tuple, T_container>::operator++()
{
    // Usual stack-style increment
    const int Nminus1 = N_containers - 1;

    int i = Nminus1;

    // Short-circuit for most common case
    // (just increment the last iterator)

    if((++iters_[i]) != enditers_[i])
    {
        tuple_[i] = *iters_[i];
        return;
    }

    // Less common cases

    for (--i; i >= 0; --i)
    {
        ++iters_[i];
        if (iters_[i] != enditers_[i])
            break;
    }

    if (i == -1)
    {
        endflag_ = true;
        return;
    }

    tuple_[i] = *iters_[i];

    for (++i; i < N_containers; ++i)  
    {
        iters_[i] = firstiters_[i];
        tuple_[i] = *iters_[i];
    }
}

}

#endif // BZ_ARRAY_CARTESIAN_H

