#ifndef BLUENOISE_H
#define BLUENOISE_H

#include <array>
typedef std::array<unsigned int, 3> VecInt;

/// util.h:
#include <vector>

template<class T>
inline T sqr(const T &x)
{ return x*x; }

// transforms even the sequence 0,1,2,3,... into reasonably good random numbers
inline unsigned int randhash(unsigned int seed){
	unsigned int i=(seed^12345391u)*2654435769u;
	i^=(i<<6)^(i>>26);
	i*=2654435769u;
	i+=(i<<5)^(i>>12);
	return i;
}

// returns repeatable stateless pseudo-random number in [0,1]
inline double randhashd(unsigned int seed){ return randhash(seed)/(double)UINT_MAX; }
inline float randhashf(unsigned int seed){ return randhash(seed)/(float)UINT_MAX; }

// returns repeatable stateless pseudo-random number in [a,b]
inline double randhashd(unsigned int seed, double a, double b){ return (b-a)*randhash(seed)/(double)UINT_MAX + a; }
inline float randhashf(unsigned int seed, float a, float b){ return ( (b-a)*randhash(seed)/(float)UINT_MAX + a); }

template<class T>
void erase_unordered(std::vector<T> &a, unsigned int index){
	a[index]=a.back();
	a.pop_back();
}
/// End of 'util.h'

/// vec.h:
template<unsigned int N, class T, class Vec>
inline T dist2(const Vec &a, const Vec &b){ 
	T d=sqr(a[0]-b[0]);
	for(unsigned int i=1; i<N; ++i) d+=sqr(a[i]-b[i]);
	return d;
}
template<unsigned int N, class T, class Vec>
T mag2(const Vec &a){
	T l=sqr(a[0]);
	for(unsigned int i=1; i<N; ++i) l+=sqr(a[i]);
	return l;
}
/// End of 'vec.h'

template<unsigned int N, class T, class Vec>
void sample_annulus(T radius, const Vec &centre, unsigned int &seed, Vec &x)
{
	Vec r;
	for(;;){
		for(unsigned int i=0; i<N; ++i){
			r[i]=4*(randhash(seed++)/(T)UINT_MAX-(T)0.5);
		}
		T r2=mag2<N,T,Vec>(r);
		if(r2>1 && r2<=4)
			break;
	}
	x=centre+radius*r;
}

template<unsigned int N, class T, class Vec>
inline Vec toVec3(VecInt j){
	Vec v;
	for(unsigned int i=0; i<N; ++i) v[i] = j[i];
	return v;
}

template<unsigned int N, class T, class Vec>
unsigned long int n_dimensional_array_index(const VecInt &dimensions, const Vec &x)
{
   unsigned long int k=0;
   if(x[N-1]>=0){
      k=(unsigned long int)x[N-1];
      if(k>=dimensions[N-1]) k=dimensions[N-1]-1;
   }
   for(unsigned int i=N-1; i>0; --i){
      k*=dimensions[i-1];
      if(x[i-1]>=0){
         unsigned int j=(int)x[i-1];
         if(j>=dimensions[i-1]) j=dimensions[i-1]-1;
         k+=j;
      }
   }
   return k;
}

template<unsigned int N, class T, class Vec>
void bluenoise_sample(T radius, Vec xmin, Vec xmax, std::vector<Vec > &sample,
                      unsigned int seed=0, int max_sample_attempts=30)
{
   sample.clear();
   std::vector<size_t> active_list;

   // acceleration grid
   T grid_dx=T(0.999)*radius/std::sqrt((T)N); // a grid cell this size can have at most one sample in it
   VecInt dimensions;
   unsigned long int total_array_size=1;
   for(unsigned int i=0; i<N; ++i){
      dimensions[i]=(unsigned int)std::ceil((xmax[i]-xmin[i])/grid_dx);
      total_array_size*=dimensions[i];
   }
   std::vector<int> accel(total_array_size, -1); // -1 indicates no sample there; otherwise index of sample point

   // first sample
   Vec x;
   for(unsigned int i=0; i<N; ++i){
      x[i]=(xmax[i]-xmin[i])*(randhash(seed++)/(T)UINT_MAX) + xmin[i];
   }
   sample.push_back(x);
   active_list.push_back(0);
   unsigned int k=n_dimensional_array_index<N,T,Vec>(dimensions, (x-xmin)/grid_dx);
   accel[k]=0;

   while(!active_list.empty()){
      unsigned int r=int(randhashf(seed++, 0, active_list.size()-0.0001f));
      size_t p=active_list[r];
      bool found_sample=false;
      VecInt j, jmin, jmax;
      for(int attempt=0; attempt<max_sample_attempts; ++attempt){
         sample_annulus<N>(radius, sample[p], seed, x);
         // check this sample is within bounds
         for(unsigned int i=0; i<N; ++i){
            if(x[i]<xmin[i] || x[i]>xmax[i])
               goto reject_sample;
         }
         // test proximity to nearby samples
         for(unsigned int i=0; i<N; ++i){
            int thismin=(int)((x[i]-radius-xmin[i])/grid_dx);
            if(thismin<0) thismin=0;
            else if(thismin>=(int)dimensions[i]) thismin=dimensions[i]-1;
            jmin[i]=(unsigned int)thismin;
            int thismax=(int)((x[i]+radius-xmin[i])/grid_dx);
            if(thismax<0) thismax=0;
            else if(thismax>=(int)dimensions[i]) thismax=dimensions[i]-1;
            jmax[i]=(unsigned int)thismax;
         }
         for(j=jmin;;){
            // check if there's a sample at j that's too close to x
            k=n_dimensional_array_index<N,T,Vec>(dimensions, toVec3<N,T,Vec>(j));
            if(accel[k]>=0 && accel[k]!=p){ // if there is a sample point different from p
               if(dist2<N,T>(x, sample[accel[k]])<sqr(radius))
                  goto reject_sample;
            }
            // move on to next j
            for(unsigned int i=0; i<N; ++i){
               ++j[i];
               if(j[i]<=jmax[i]){
                  break;
               }else{
                  if(i==N-1) goto done_j_loop;
                  else j[i]=jmin[i]; // and try incrementing the next dimension along
               }
            }
         }
         done_j_loop:
         // if we made it here, we're good!
         found_sample=true;
         break;
         // if we goto here, x is too close to an existing sample
         reject_sample:
         ; // nothing to do except go to the next iteration in this loop
      }
      if(found_sample){
         size_t q=sample.size(); // the index of the new sample
         sample.push_back(x);
         active_list.push_back(q);
         k=n_dimensional_array_index<N,T,Vec>(dimensions, (x-xmin)/grid_dx);
         accel[k]=(int)q;
      }else{
         // since we couldn't find a sample on p's disk, we remove p from the active list
         erase_unordered(active_list, r);
      }
   }
}

#endif
