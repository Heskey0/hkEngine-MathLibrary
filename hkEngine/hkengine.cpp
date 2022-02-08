
#include "math.hpp"



using namespace hkEngine;

int main(int argc, char* argv[]) {
	float3 c0 = float3(1, 1, 1);
	Transform t1 = Translate(float3(2, 2, 2));
	float3 c1 = t1(c0);
	
	std::cout << t1 << std::endl;
	std::cout << c0 << c1 <<std::endl;
}