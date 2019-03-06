#pragma once
#include <ostream>
using namespace std;
/*
如果不加名字空间，直接写operator<<的重载，编译器只会尝试匹配std::operator<<，找不到则提示error，如果这样写
template <typename U>
friend ostream & operator<<(ostream &out, const Widget &a);
其中Widget是一个类，而不是一个类模板，则编译器会提示为模板参数U推导类型
另外所有使用ARRAY::operator << 的地方都要显示的写ARRAY::operator << 的定义
*/
namespace ARRAY
{
	template <typename T>
	class Array
	{
	public:
		template <typename T>
		friend ostream & operator<<(ostream &out, const T &a);
	private:
	};
}

