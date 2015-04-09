/**
* Code by: Robin Green and others
*/

#pragma once
#include <vector>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265f
#endif

template<typename Vec3d, typename Scalar>
struct SHSample 
{ 
	SHSample( int numBands = 4 ){ sph = vec = Vec3d(0,0,0); coeff.resize(numBands*numBands,0); }
	Vec3d sph;
	Vec3d vec;
	std::vector<Scalar> coeff; 
}; 

template<typename Vec3d, typename Scalar>
class SphericalHarmonic
{
public:
	SphericalHarmonic( int numBands = 4 ){ 
		generateFactorials( 20 ); 
		init( numBands );
	}

	void init( int numBands ){ 
		n_bands = numBands; 
		n_coeff = numBands*numBands; 
	}

	long long generateFactorial(int n){
		if(n > 0)
			return n * generateFactorial(n - 1);
		else
			return 1;
	}

	void generateFactorials(int n){
		for(int i = 0; i <= n; ++i){
			if( i > 0){
				long long v = i * _factorial[i - 1];
				_factorial.push_back(v);
			}
			else{
				_factorial.push_back(1);
			}
		}
	}

	Scalar K(int l, int m) 
	{ 
		// re-normalization constant for SH function 
		Scalar temp = ((2.0*l+1.0) * _factorial[l-m]) / (4.0 * M_PI * _factorial[l+m]); 
		return sqrt(temp); 
	}

	Scalar P(int l,int m, Scalar x) 
	{ 
		// evaluate an Associated Legendre Polynomial P(l,m,x) at x 
		Scalar pmm = 1.0; 
		if(m>0) { 
			Scalar somx2 = sqrt((1.0-x)*(1.0+x)); 
			Scalar fact = 1.0; 
			for(int i=1; i<=m; i++) { 
				pmm *= (-fact) * somx2; 
				fact += 2.0; 
			} 
		} 
		if(l==m) 
		{
			return pmm; 
		}
		Scalar pmmp1 = x * (2.0*m+1.0) * pmm; 
		if(l==m+1) 
		{
			return pmmp1; 
		}
		Scalar pll = 0.0; 
		for(int ll=m+2; ll<=l; ++ll) { 
			pll = ( (2.0*ll-1.0)*x*pmmp1-(ll+m-1.0)*pmm ) / (ll-m); 
			pmm = pmmp1; 
			pmmp1 = pll; 
		} 
		return pll; 
	}

	Scalar SH(int l, int m, Scalar theta, Scalar phi) 
	{ 
		// return a point sample of a Spherical Harmonic basis function 
		// l is the band, range [0..N] 
		// m in the range [-l..l] 
		// theta in the range [0..Pi] 
		// phi in the range [0..2*Pi] 
		const Scalar sqrt2 = sqrt(2.0); 
		if(m == 0) 
			return K(l,0) * P(l,m,cos(theta)); 
		else if(m > 0) 
			return sqrt2*K(l,m) * cos(m*phi) * P(l,m,cos(theta)); 
		else 
			return sqrt2*K(l,-m) * sin(-m*phi) * P(l,-m,cos(theta)); 
	}

	void SH_setup_spherical_samples_uniform(std::vector< SHSample<Vec3d, Scalar> >& samples, int sqrt_n_samples) 
	{ 
		samples.resize(sqrt_n_samples * sqrt_n_samples);
		// fill an N*N*2 array with uniformly distributed 
		// samples across the sphere using jitter-ed stratification 
		int i=0;  // array index 
		Scalar oneoverN = 1.0/sqrt_n_samples; 
		for(int a=0; a < sqrt_n_samples; a++) { 
			for(int b=0; b < sqrt_n_samples; b++) { 
				// generate unbiased distribution of spherical coordinates
				samples[i] = SHSample<Vec3d>( n_bands );
				Scalar x = (a + rand()/(Scalar)RAND_MAX) * oneoverN;  // do not reuse results 
				Scalar y = (b + rand()/(Scalar)RAND_MAX) * oneoverN;  // each sample must be random 
				Scalar theta = 2.0 * acos(sqrt(1.0 - x)); 
				Scalar phi = 2.0 * M_PI * y;
				samples[i].sph = Vec3d(theta,phi,1.0); 
				// convert spherical coordinates to unit vector 
				Vec3d vec(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)); 
				samples[i].vec = vec; 
				// pre-compute all SH coefficients for this sample 
				for(int l=0; l < n_bands; ++l) { 
					for(int m=-l; m<=l; ++m) { 
						int index = l*(l+1)+m; 
						samples[i].coeff[index] = SH(l,m,theta,phi); 
					} 
				} 
				++i; 
			} 
		} 
	} 

	void SH_setup_spherical( const std::vector<Vec3d> & directions, std::vector< SHSample<Vec3d,Scalar> >& samples )
	{
		samples.resize( directions.size() );
		int i = 0;

		for(auto d : directions)
		{
			samples[i] = SHSample<Vec3d,Scalar>( n_bands );
			samples[i].vec = d; 

			Scalar r = std::sqrt(pow(d.x(),2) + pow(d.y(),2) + pow(d.z(),2));
			Scalar theta = acos( d.z() / r );
			Scalar phi = atan2(d.y(), d.x());
			samples[i].sph = Vec3d(theta,phi,1.0); 

			// pre-compute all SH coefficients for this sample 
			for(int l=0; l < n_bands; ++l) { 
				for(int m=-l; m<=l; ++m) { 
					int index = l*(l+1)+m; 
					samples[i].coeff[index] = SH(l,m,theta,phi); 
				}
			}

			i++;
		}
	}

	typedef Scalar (*SH_polar_fn)(Scalar theta, Scalar phi); 

	void SH_project_polar_function(SH_polar_fn fn, const std::vector< SHSample<Vec3d,Scalar> >& samples, std::vector<Scalar>& result) 
	{ 
		result.clear();
		result.resize(n_coeff,0);

		const Scalar weight = 4.0 * M_PI; 
		int n_samples = samples.size();

		// for each sample 
		for(int i=0; i < n_samples; ++i) 
		{ 
			Scalar theta = samples[i].sph.x(); 
			Scalar phi   = samples[i].sph.y(); 
			for(int n=0; n < n_coeff; ++n) 
			{ 
				result[n] += fn(theta, phi) * samples[i].coeff[n]; 
			}
		} 

		// divide the result by weight and number of samples 
		Scalar factor = weight / n_samples; 
		for(int i=0; i < n_coeff; ++i) 
		{ 
			result[i] = result[i] * factor; 
		} 
	}

	void SH_project_function(const std::vector<Scalar> & fn_values, const std::vector< SHSample<Vec3d,Scalar> >& samples, std::vector<Scalar>& result) const
	{ 
		result.clear();
		result.resize(n_coeff, 0);

		const Scalar weight = Scalar(4.0 * M_PI); 
		size_t n_samples = samples.size();

		// for each sample 
		for(size_t i=0; i < n_samples; ++i) 
			for(size_t n=0; n < n_coeff; ++n)
				result[n] += fn_values[i] * samples[i].coeff[n]; 

		// divide the result by weight and number of samples 
		Scalar factor = weight / n_samples; 
		for(size_t i=0; i < n_coeff; ++i)
			result[i] = result[i] * factor;
	}

	inline std::vector<Scalar> SH_signature( const std::vector<Scalar> & coeffs  ) const
	{
		std::vector<Scalar> output( n_bands, 0.0 );

		for(int l=0; l < n_bands; ++l) 
		{
			Scalar norm2 = 0;

			for(int m = -l; m <= l; ++m)
			{
				int index = l*(l+1)+m; 
				Scalar r = coeffs[index];

				if(index == 0)
					norm2 += (r*r);
				else
					norm2 += (r*r)*2;
			}

			output[l] = std::sqrt( norm2 );
		}

		return output;
	}

	std::vector<Scalar> SH_reconstruct( const std::vector<Vec3d,Scalar> & positions, std::vector<Scalar> & result  )
	{
		std::vector<Scalar> output;

		for(auto p: positions)
		{
			Scalar r = std::sqrt(pow(p.x(),2) + pow(p.y(),2) + pow(p.z(),2));
			Scalar theta = acos( p.z() / r );
			Scalar phi = atan2(p.y(), p.x());

			Scalar val = 0;
			for(int l=0; l<n_bands; ++l) {
				for(int m=-l; m<=l; ++m) {
					int index = l*(l+1)+m;
					val += result[index] * SH(l,m,theta,phi);
				}
			}
			output.push_back( std::max(0.0, val) );
		}

		return output;
	}

private:
	int n_bands;
	int n_coeff;
	std::vector<long long> _factorial;
};

