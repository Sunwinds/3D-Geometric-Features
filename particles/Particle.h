#pragma once

#include <vector>
#include <stdint.h>
#include <Eigen/Core>

class Particle
{
public:
    Particle(const Eigen::Vector3d& pos) : pos(pos), measure(0.0), weight(1.0), alpha(1.0), direction(Eigen::Vector3d(0,0,1)), flag(0) {}

    size_t id, correspondence;
	int flag;
	uint64_t morton;
    Eigen::Vector3d pos, direction, relativePos;
    double measure;
	double weight;
	double alpha;
};
