#pragma once
#include <ostream>
using namespace std;
/*
����������ֿռ䣬ֱ��дoperator<<�����أ�������ֻ�᳢��ƥ��std::operator<<���Ҳ�������ʾerror���������д
template <typename U>
friend ostream & operator<<(ostream &out, const Widget &a);
����Widget��һ���࣬������һ����ģ�壬�����������ʾΪģ�����U�Ƶ�����
��������ʹ��ARRAY::operator << �ĵط���Ҫ��ʾ��дARRAY::operator << �Ķ���
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

