#pragma once
#include <iostream>

#define HK_CHECK(condition) if(!condition)\
{\
	std::cout << std::endl << std::endl;\
	std::cout << "!! Error::" #condition << std::endl;\
	std::cout << "File::" << __FILE__ << std::endl;\
	std::cout << "Line::" << __LINE__ << std::endl;\
	Error(#condition);\
	std::cout << std::endl << std::endl;\
}

inline void Error(std::string info) {
	std::cout << "[#############################]" << std::endl;
}
