#include "hiradar/interface.hpp" // 包含对应的头文件

// 为 IPowerDensityWriter 的虚析构函数提供一个空的实现体
IPowerDensityWriter::~IPowerDensityWriter() {}

// 为 IRadar 的虚析构函数提供一个空的实现体
// 这个函数体虽然是空的，但它为编译器生成vtable提供了“锚点”，从而解决链接错误
IRadar::~IRadar() {}