#ifndef parameters_hh
#define parameters_hh

namespace param
{
	const float gravity = 9.81f;  //direction of the gravity is in the positive y-axis (down)
	const float dt = 0.1f;         //time step in the calculation of shots etc. increasing dt makes bullets go faster
    const float windEffect = 1.0f; // coefficient for effect of wind on projectile
}

#endif
