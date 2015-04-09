#pragma once

// area3D_Polygon(): compute the area of a 3D planar polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+1 points in a 2D plane with V[n]=V[0]
//          Point N = a normal vector of the polygon's plane
//  Return: the (double) area of the polygon
template <class Point>
double area3D_Polygon( const std::vector<Point> & V, Point * normal = NULL )
{
    double area = 0;
    double an, ax, ay, az; // abs value of normal and its coords
    int  coord;           // coord to ignore: 1=x, 2=y, 3=z
    int  i, j, k;         // loop indices

	int n = V.size();

    if (n < 3) return 0;  // a degenerate polygon

	Point N(0,0,1);
	if(!normal)
	{	
		// Compute normal if not provided
		Point mean(0,0,0);
		for(auto p : V) mean += p;
		mean /= V.size();
		N = (V[0] - mean).cross(V[1] - mean).normalized();
	}
	else
	{
		N = *normal;
	}

    // select largest abs coordinate to ignore for projection
    ax = (N.x() > 0 ? N.x() : -N.x());    // abs x-coord
    ay = (N.y() > 0 ? N.y() : -N.y());    // abs y-coord
    az = (N.z() > 0 ? N.z() : -N.z());    // abs z-coord

    coord = 3;                    // ignore z-coord
    if (ax > ay) {
        if (ax > az) coord = 1;   // ignore x-coord
    }
    else if (ay > az) coord = 2;  // ignore y-coord

    // compute area of the 2D projection
    switch (coord) {
      case 1:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].y() * (V[j].z() - V[k].z()));
        break;
      case 2:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].z() * (V[j].x() - V[k].x()));
        break;
      case 3:
        for (i=1, j=2, k=0; i<n; i++, j++, k++)
            area += (V[i].x() * (V[j].y() - V[k].y()));
        break;
    }
    switch (coord) {    // wrap-around term
      case 1:
        area += (V[n].y() * (V[1].z() - V[n-1].z()));
        break;
      case 2:
        area += (V[n].z() * (V[1].x() - V[n-1].x()));
        break;
      case 3:
        area += (V[n].x() * (V[1].y() - V[n-1].y()));
        break;
    }

    // scale to get area before projection
    an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
    switch (coord) {
      case 1:
        area *= (an / (2 * N.x()));
        break;
      case 2:
        area *= (an / (2 * N.y()));
        break;
      case 3:
        area *= (an / (2 * N.z()));
    }
    return area;
}
//===================================================================
