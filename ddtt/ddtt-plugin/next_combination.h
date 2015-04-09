// Adapted from code by Fredrik Eriksson
#pragma once
#include <algorithm>

namespace std
{
	template<typename BidirectionalIterator>
	bool next_combination(BidirectionalIterator first, BidirectionalIterator k, BidirectionalIterator last){
		if( first == last || k == first || k == last){
			return false;
		}
		BidirectionalIterator i = first;
		BidirectionalIterator ii = last;

		++i;
		if( i == last ){
			return false;
		}

		--i;
		i = k;
		--ii;

		while ( i != first ){
			if ( *--i < *ii ){
				BidirectionalIterator j = k;
				while ( !(*i < *j) ){
					++j;
				}
				std::iter_swap( i, j );
				++i;
				++j;
				ii = k;
				std::rotate( i, j, last );
				while ( j != last ){
					++j;
					++ii;
				}
				std::rotate( k, ii, last );
				return true;
			}
		}
		std::rotate( first, k, last );
		return false;
	}

	template<class InputContainer, class OutputContainer>
	size_t all_combinations(InputContainer input, int k, OutputContainer & output){
		InputContainer::iterator first = input.begin(), last = input.end();
		do{
			InputContainer comb(first, first + k);
			output.push_back(comb);
		} while( std::next_combination(first, first + k, last ) );
		return output.size();
	}
}
