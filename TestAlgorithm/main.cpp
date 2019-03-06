#include <stdlib.h>
#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <list>
#include <map>
#include <time.h>
#include <functional>
#include <random>
#include <limits>
using namespace std;

#include "Array.h"

void print(int *a, int count)
{
    for(int i = 0; i < count; ++i)
        cout << a[i] << " ";
    cout << endl;
}

template<typename T>
void swap(T &a, T &b)
{
    T temp = a;
    a = b;
    b = temp;
}

// simple_search
int simple_search(int *a, int count, int key)
{
    for(int i = 0; i < count; ++i)
    {
        if(key == a[i])
            return i;
    }
	return -1;
}

// binary_search
int binary_search(int *a, int count, int key)
{
    int low = 0, high = count - 1;
    int mid = 0;
    while(low <= high)
    {
        mid = low + (high - low) / 2;
        if(key == a[mid])
            return mid;
		else if (key > a[mid])
            low = mid + 1;
        else
            high = mid - 1;
    }
	return low;
}

// recursive_binary_search
int recursive_binary_search(int *a, int low, int high, int key)
{
    if(low > high) return low;
    int mid = low + (high - low) / 2;
    if(key == a[mid]) return mid;
    else if(key > a[mid]) return recursive_binary_search(a, mid+1, high, key);
    else return recursive_binary_search(a, low, mid-1, key);
}

//O(n ^ 2)
void simple_select(int *a, int count)
{
    for(int i = 0; i < count - 1; ++i)
    {
        int k = i;
        for(int j = i+1; j < count; ++j)
        {
            if(a[j] < a[k])
                k = j;
        }
		::swap(a[i], a[k]);
    }
}

// O(n^2)
void bubble_sort(int *a, int count)
{
    bool bChange = true;
    for(int i = 0; i < count - 1 && bChange; ++i)
    {
        bChange = false;
        for(int j = count - 1; j > i; --j)
        {
            if(a[j-1] > a[j])
            {
                bChange = true;
                ::swap(a[j-1], a[j]);
            }
        }
    }
}

// O(n^2)
void insert_sort(int *a, int count)
{
    for(int i = 1; i < count; ++i)
    {
        int key = a[i];
        int j = i - 1;
        while(j >= 0 && a[j] > key)
        {
            a[j+1] = a[j];
            --j;
        }
        a[j+1] = key;
    }
}

// quick_sort O(nlgn)
int do_quick_sort(int *a, int p, int q)
{
    int s = p, t = q;
    int key = a[s];
    while(s < t)
    {
        while(s < t && a[t] >= key) --t;
        a[s] = a[t];
        while(s < t && a[s] <= key) ++s;
        a[t] = a[s];
    }
    a[s] = key;
    return s;
}

int _do_quick_sort(int *a, int p, int q)
{
	int key = a[p];
	int s = p - 1, t = q + 1;
	while (true)
	{
		do 
		{
			--t;
		} while (a[t] > key);
		do
		{
			++s;
		} while (a[s] < key);
		if (s < t)
			::swap(a[s], a[t]);
		else return s;
	}
}

int __do_quick_sort(int *a, int p, int q)
{
	int key = a[q];
	int i = p - 1;
	for (int j = p; j < q; ++j)
	{
		if (a[j] < key)
		{
			++i;
			::swap(a[i], a[j]);
		}
	}
	::swap(a[++i], a[q]);
	return i;
}

void quick_sort(int *a, int p, int q)
{
    if(p >= q) return;
	int k = _do_quick_sort(a, p, q);
    quick_sort(a, p, k - 1);
    quick_sort(a, k + 1, q);
}

//merge_sort O(nlgn)
void merge_array(int *a, int p, int k, int q)
{
	//merge a[p...k] and a[k+1...q] to a[p...q]
    int n1 = k - p + 1;
    int n2 = q - (k + 1) + 1;
    int *a1 = new int[n1];
    int *a2 = new int[n2];
    memcpy(a1, a + p, n1 * sizeof(int));
	memcpy(a2, a + k + 1, n2 * sizeof(int));
    int i = 0, j = 0, n = p;
    while(i < n1 && j < n2)
    {
        if(a1[i] <= a2[j])
            a[n++] = a1[i++];
        else
            a[n++] = a2[j++];
    }
    if(i >= n1)
    {
        while(j < n2)
            a[n++] = a2[j++];
    }

    if(j >= n2)
    {
        while(i < n1)
            a[n++] = a1[i++];
    }
    delete[] a1;
    delete[] a2;
}

void merge_sort(int *a, int p, int q)
{
    if(p >= q) return;
    int k = (p + q) / 2;
    merge_sort(a, p, k);
    merge_sort(a, k+1, q);
    merge_array(a, p, k, q);
}

//binary_search_insert_sort O(nlgn)
void binary_search_insert_sort(int *a, int count)
{
    for(int i = 1; i < count; ++i)
    {
        int key = a[i];
        int j = i - 1;
		int k = binary_search(a, i, key);
        while(j >= k)
        {
            a[j+1] = a[j];
            --j;
        }
        a[k] = key;
    }
}

struct subarray_result
{
	int low;
	int hight;
	int sum;
};

subarray_result get_max_subarray_cross_mid(const int *a, int low, int mid, int high)
{
	int left_index = mid, left_current_sum = 0, left_sum = 0;
	int right_index = mid + 1, right_current_sum = 0, right_sum = 0;
	for (int i = mid; i >= low; --i)
	{
		left_current_sum += a[i];
		if (left_current_sum > left_sum)
		{
			left_sum = left_current_sum;
			left_index = i;
		}
	}

	for (int j = mid + 1; j <= high; ++j)
	{
		right_current_sum += a[j];
		if (right_current_sum > right_sum)
		{
			right_sum = right_current_sum;
			right_index = j;
		}
	}
	return { left_index, right_index, left_sum + right_sum };
}
subarray_result get_max_subarray(const int *a, int low, int high)
{
	subarray_result result_left{ low, high, a[low] };
	subarray_result result_right{ low, high, a[low] };
	subarray_result result_cross{ low, high, a[low] };
	if (low == high)
		return subarray_result{ low, high, a[low] };
	int mid = (low + high) / 2;
	result_left = get_max_subarray(a, low, mid);
	result_right = get_max_subarray(a, mid + 1, high);
	result_cross = get_max_subarray_cross_mid(a, low, mid, high);
	if (result_left.sum >= result_right.sum && result_left.sum >= result_cross.sum)
		return result_left;
	if (result_right.sum >= result_left.sum && result_right.sum >= result_cross.sum)
		return result_right;
	return result_cross;
}

struct tHeap
{
	int parent(int i){ return i / 2; }
	int left(int i){ return 2 * i + 1; }
	int right(int i){ return 2 * i + 2; }
	void max_heap(int i)
	{
		int l = left(i);
		int r = right(i);
		int largest = i;
		if (l < heap_size && a[l] > a[largest])
			largest = l;
		if (r < heap_size && a[r] > a[largest])
			largest = r;
		if (i != largest)
		{
			::swap(a[i], a[largest]);
			max_heap(largest);
		}
	}
	void build_heap()
	{
		heap_size = length;
		for (int i = (length - 1) / 2; i >= 0; --i)
			max_heap(i);
	}

	void heap_sort()
	{
		for (int i = length - 1; i >= 0; --i)
		{
			::swap(a[0], a[i]);
			--heap_size;
			max_heap(0);
		}
	}

	int length;
	int heap_size;
	int a[1024];
};

void count_sort(int *a, int length, int k)
{
	int *b = new int[length];
	memcpy(b, a, length * sizeof(int));
	int *c = new int[k];
	memset(c, 0, k * sizeof(int));
	memset(a, 0, length * sizeof(int));
	for (int i = 0; i < length; ++i)
		++c[b[i]];

	for (int i = 1; i < k; ++i)
		c[i] += c[i - 1];

	for (int i = length - 1; i >= 0; --i)
	{
		a[c[b[i]]-1] = b[i];
		--c[b[i]];
	}

	delete[] c;
	delete[] b;
}

//random

double random_0_1()
{
	random_device rd;
	default_random_engine gen(rd());
	uniform_real_distribution<> dis(0, 1);
	return dis(gen);
}

int random_int()
{
	random_device rd;
	default_random_engine gen(rd());
	uniform_int_distribution<> dis(0, std::numeric_limits<int>::max());
	return dis(gen);
}

double uniform_real(double a, double b)
{
	return a + random_0_1() * (b - a);
}

int uniform_int(int a, int b)
{
	return static_cast<int>(uniform_real(a, b));
}

int uniform_n(int n)
{
	return uniform_int(0, n);
}

void shuffle(int *a, int n)
{
	for (int i = 0; i < n; ++i)
	{
		std::swap(a[i], a[uniform_int(i, n)]);
	}
}

void randomSelect(const int *src, int *dst, int k, int n)
{
	int *src_t = new int[n];
	memcpy(src_t, src, n * sizeof(int));
	/*
	auto i = 0, temp = 0;
	for (int j = 0; j < k; ++j)
	{
		i = random_int() % n;
		temp = src_t[i];
		src_t[i] = src_t[j];
		src_t[j] = temp;
	}
	i = 0;
	*/
	shuffle(src_t, n);
	for (int i = 0; i < k; ++i)
		dst[i] = src_t[i];
	delete[] src_t;
}

// p(i) = a[i]
// i = 0, p(0) => 0 <= r <= a[0] => a[0] - 0 = a[0]
//if p(i) => S(0...i-1) <= r <= S(0...i) => S(0...i) - S(0...i-1) = a[i]
void discrete(double *a, int n)
{
	double r = random_0_1();
	double sum = 0;
	for (int i = 0; i < n; ++i)
	{
		sum += a[i];
		if (sum >= r) return i;
	}
	return -1;
}

// dp
// 20 floors, only 1 step or 2 steps every time, what is the max total step method?
// from up to down, memory
// the question is separated to 2 questions with two choices, a: the last step is one step, then the sub question is what the max total step is the rest 19 floors.
// b: the last step is two step, then the sub question is what the max total steps is the rest 18 floors.
//  n > 0
//  S(n) = n, when n <= 2
//  S(n) = S(n-1) + S(n-2)
int max_step(int *a, int n)
{
	if (a[n] <= 0)
		a[n] = max_step(a, n - 1) + max_step(a, n - 2);
	return a[n];
}

int memory_max_step(int n)
{
	int *a = new int[n];
	memset(a, 0, n * sizeof(int));
	a[0] = 0;
	a[1] = 1;
	a[2] = 2;
	int ret = max_step(a, n);
	return ret;
}

// there are 9 candys, at least 1 will be eaten every day, what is the total eating method?
// the last day eat 1, 2, 3, ... , 8, 9,then the sub question is what is the total eating method if there are 8, 7, 6, ..., 1, 0
//  n >= 0
//  S(n) = 1, when n <= 1
//  S(n) = S(n-1) + S(n-2) + ... + S(n-n)
int _max_eat(int *a, int n)
{
	if (a[n] <= 0)
	{
		int result = 0;
		for (int i = 1; i <= n; ++i)
		{
			result += _max_eat(a, n - i);
		}
		a[n] = result;
	}
	return a[n];
}

int _memory_max_eat(int n)
{
	int *a = new int[n + 1];
	memset(a, 0, (n + 1) * sizeof(int));
	a[0] = 1;
	a[1] = 1;
	int ret = _max_eat(a, n);
	delete[] a;
	return ret;
}
// what if the candy is different?
// 同样用之前的那个256种分隔方案，假设某天要吃n颗糖，当前剩余m颗糖，那么就从剩余的m颗糖中选择n颗，有C(n, m)种选法，伪代码是
// sum = 0
// for i = 1, 256 do
// 假设e[i]为第i个分隔方案，e[i][j]为第i个分隔方案的第j天要吃的糖的数量，left(e[i], j)求得e[i]的第j天还剩多少糖
// result = 1
// for j = 1, e[i].size do
// result *= C(e[i][j]，left(e[i], j))
// end
// sum += result
// end

int cmn(int m, int n)
{
	int a = 1;
	int b = 1;
	for (int i = m; i >= m - n + 1; --i)
	{
		a *= i;
	}
	for (int i = n; i >= 1; --i)
	{
		b *= i;
	}
	return a / b;
}

int max_eat(int *a, int n)
{
	if (a[n] <= 0)
	{
		int result = 0;
		for (int i = 1; i <= n; ++i)
		{
			result += cmn(n, i) * max_eat(a, n - i);
		}
		a[n] = result;
	}
	return a[n];
}

int memory_max_eat(int n)
{
	int *a = new int[n + 1];
	memset(a, 0, (n + 1) * sizeof(int));
	a[0] = 1;
	a[1] = 1;
	int ret = max_eat(a, n);
	delete[] a;
	return ret;
}

/*
function cnr(n, r, s, t)
	if #t == r then
		local result = ""
		for _, v in ipairs(t) do
			result = result .. v .. " "
		end
		print(result .. ",")
		return
	end
	for i=s, #n do
		table.insert(t, n[i])
		cnr(n, r, i+1, t)
		table.remove(t, #t)
	end
end

function swap(n, i, j)
	if i == j then return end
	local tmp = n[i]
	n[i] = n[j]
	n[j] = tmp
end

function anr(n, r, s, e, t)
	--[[
	if #t == r then
		local result = ""
		for _, v in ipairs(t) do
			result = result .. v .. " "
		end
		print(result .. ",")
		return
	end
	]]
	if s >= e then
		printT(n)
		return
	end

	for i=s, e do
		swap(n, i, s)
		--table.insert(t, n[i])
		anr(n, r, i+1, e, t)
		--table.remove(t, #t)
		swap(n, i, s)
	end
end

local count = 0
function sumof(n, c, s, t)
	if n < 0 then return end
		if 0 == n and #t == c then
		local result = ""
		for _, v in ipairs(t) do
			result = result .. v .. " "
		end
		print(result .. ",")
		count = count + 1
		return
	end

	for i=s, n do
		table.insert(t, i)
		sumof(n-i, c, i+1, t)
		table.remove(t, #t)
	end
end

local t = {}
sumof(10, 4, 1, t)
*/

void SumOf(int n, int s, std::vector<int> &t)
{
	if (0 == n)
	{
		for (auto v : t)
		{
			printf("%d ", v);
		}
		cout << endl;
	}

	for (int i = s; i <= n; ++i)
	{
		t.push_back(i);
		SumOf(n - i, i+1, t);
		t.pop_back();
	}
}

void Cmn(int *n, int r, int s, int e, std::vector<int> &t)
{
	if (r == t.size())
	{
		for (auto v : t)
		{
			printf("%d ", v);
		}
		cout << endl;
	}

	for (int i = s; i <= e; ++i)
	{
		t.push_back(n[i]);
		Cmn(n, r, i + 1, e, t);
		t.pop_back();
	}
}

void Amn(int *n, int s, int e)
{
	if (s == e)
	{
		for (int i = 0; i <= e; ++i)
		{
			printf("%d", n[i]);
		}
		cout << endl;
		return;
	}

	for (int i = s; i <= e; ++i)
	{
		std::swap(n[s], n[i]);
		Amn(n, s + 1, e);
		std::swap(n[s], n[i]);
	}
}

void Permutation(char* pStr, char* pBegin)
{
	if (!pStr || !pBegin)
		return;

	if (*pBegin == '\0')
	{
		printf("%s\n", pStr);
	}
	else
	{
		for (char* pCh = pBegin; *pCh != '\0'; ++pCh)
		{
			// swap pCh and pBegin    
			char temp = *pCh;
			*pCh = *pBegin;
			*pBegin = temp;

			Permutation(pStr, pBegin + 1);
			// restore pCh and pBegin    
			temp = *pCh;
			*pCh = *pBegin;
			*pBegin = temp;
		}
	}
}


bool compare_without_anything(int *a, int length)
{
	int *p = a, *q = a + length - 1;
	return (p < q ? 
		p[1] < p[0] ? false : compare_without_anything(a + 1, length - 1) 
		: true);
}

unsigned long fibonacci_n()
{
	static unsigned long a = 0L;
	static unsigned long b = 1L;
	unsigned long c = a + b;
	a = b;
	b = c;
	return a;
}

namespace ARRAY
{
	template <typename T>
	ostream & operator<<(ostream &out, const T &a)
	{
		return out;
	}
}

//Graghic
const int n = 6;
struct Vertex
{
	char c;
	int idx;
	int color;//white:0, gray:1, black:2
	int d;
};

void InitGraphic(Vertex *g, int s)
{
	for (auto i = 0; i < n; ++i)
	{
		g[i].color = 0;
		g[i].d = INT_MAX;
	}
	g[s].color = 1;
	g[s].d = 0;
}

void BFS(Vertex *g, int (*gm)[n], int s)
{
	InitGraphic(g, s);
	std::queue<Vertex> q;
	q.push(g[s]);
	cout << "The BFS of the Graphic is:" << endl;
	while (!q.empty())
	{
		auto &v = g[q.front().idx];
		q.pop();
		printf("%c ", v.c);
		for (auto i = 0; i < n; ++i)
		{
			auto &u = g[i];
			if (1 == gm[v.idx][i] && 0 == u.color)
			{
				u.color = 1;
				u.d = v.d + 1;
				q.push(u);
			}
		}
		v.color = 2;
	}
}

void _DFS(Vertex *g, int(*gm)[n])
{
	g->color = 1;
	printf("%c ", g->c);
	for (auto i = 0; i < n; ++i)
	{
		auto &u = g[i];
		if (1 == gm[g->idx][i] && 0 == u.color)
		{
			_DFS(&u, gm);
		}
	}
	g->color = 2;
}

void DFS(Vertex *g, int(*gm)[n], int s)
{
	InitGraphic(g, s);
	g[s].color = 0;
	cout << endl;
	cout << "The DFS of the Graphic is:" << endl;
	for (auto i = 0; i < n; ++i)
		if (0 == g[i].color)
			_DFS(&g[i], gm);
}

namespace Graphic
{
	const int N = 5;

	struct V
	{
		enum class COLOR : int
		{
			WHITE,
			GRAY,
			BLACK,
		};

		V()
		: c('\0')
		, d(INT_MAX)
		, pi1('\0')
		, b(0)
		, f(0)
		, pi2('\0')
		, color(COLOR::WHITE)
		{

		}

		char c;
		char pi1;
		char pi2;
		char pi3;
		int d;
		int b;
		int f;
		COLOR color;
	};

	struct G1
	{
		struct V1
		{
			V v;
			std::list<V1 *> l;
		};

		G1()
		:g(N)
		{
			for (auto i = 0U; i < g.size(); ++i)
			{
				g[i].v.c = '0' + i;
			}

			g[0].l = { &g[1], &g[2], };
			g[1].l = { &g[2], &g[3], };
			g[2].l = { &g[1], &g[3], &g[4], };
			g[3].l = { &g[4], };
			g[4].l = { &g[1], &g[3], };
		}
		
		void BFS(int s)
		{
			for (auto &v : g)
			{
				v.v.color = V::COLOR::WHITE;
			}

			g[s].v.d = 0;
			g[s].v.color = V::COLOR::GRAY;

			std::queue<V1 *> q;
			q.push(&g[s]);

			while (!q.empty())
			{
				auto v = q.front();
				q.pop();
				printf("%c ", v->v.c);
				for (const auto &u : v->l)
				{
					if (u->v.color == V::COLOR::WHITE)
					{
						u->v.color = V::COLOR::GRAY;
						u->v.pi1 = v->v.c;
						u->v.d = v->v.d + 1;
						q.push(u);
					}
				}
				v->v.color = V::COLOR::BLACK;
			}
		}

		void DFS_VISIT(V1 &v, unsigned &time)
		{	
			printf("%c ", v.v.c);
			v.v.b = ++time;
			v.v.color = V::COLOR::GRAY;
			for (auto &u : v.l)
			{
				if (u->v.color == V::COLOR::WHITE)
				{
					u->v.pi2 = v.v.c;
					DFS_VISIT(*u, time);
				}
			}
			v.v.f = ++time;
			v.v.color = V::COLOR::BLACK;
		}

		void DFS()
		{
			for (auto &v : g)
			{
				v.v.color = V::COLOR::WHITE;
				v.v.b = 0;
				v.v.f = 0;
			}

			auto time = 0U;

			for (auto &v : g)
				if (v.v.color == V::COLOR::WHITE)
					DFS_VISIT(v, time);
		}

		void Dijkstra(int s)
		{
			std::function<void (int)> init = [this](int s){
				for (auto &v : g)
				{
					v.v.d = INT_MAX;
					v.v.pi3 = '\0';
				}
				g[s].v.d = 0;
				g[s].v.pi3 = '\0';
			};

			std::function<int (const V1 &, const V1 &)> w = [this](const V1 &, const V1 &)->int{
				return 0;
			};

			std::function<void (V1 &, V1 &)> relax = [this, w](V1 &u, V1 &v){
				if (v.v.d > u.v.d + w(u, v))
				{
					v.v.d = u.v.d + w(u, v);
				}
				v.v.pi3 = u.v.c;
			};

			std::vector<V1 *> S;
			auto compare = [](V1 *a, V1 *b)->bool{ return a->v.d < b->v.d; };
			std::priority_queue<V1 *, std::vector<V1 *>, decltype(compare)> Q;
			for (auto &v : g)
				Q.push(&v);

			while (!Q.empty())
			{
				auto u = Q.top();
				Q.pop();
				S.push_back(u);
				for (auto v : u->l)
					relax(*u, *v);
			}
		}

		std::vector<V1> g;
	};
}

int main()
{
	int a[] = {51, 49, 30, 45, 51, 76, 10, 3, 11, 15};
    int count = sizeof(a) / sizeof(a[0]);
    cout << "print before sort" << endl;
    print(a, count);
    //bubble_sort(a, count);
    //insert_sort(a, count);
    //quick_sort(a, 0, count-1);
    merge_sort(a, 0, count-1);
    //binary_search_insert_sort(a, count);
	//count_sort(a, count, 80);
    cout << "print after sort" << endl;
    print(a, count);

	//////////////////////////////////////////////////////////////////////////
	//2,4,-7,5,2,-1,2,-4,3}的最大子数组为{5,2,-1,2}，最大子数组的和为5+2-1+2=8
	int b[] = { 2, 4, -7, 5, 2, -1, 2, -4, 3 , 1};
	int *c = new int[10]{0};
	subarray_result max_subarray = get_max_subarray(b, 0, 9);

	tHeap heap{ 10, 10, { 5, 3, 12, 5, 2, 1, -1, 1, 2, 5 } };
	heap.build_heap();
	tHeap tmp = heap;
	tmp.heap_sort();

	int d[] = { 3,1, 2, 3,4, 6, 7, 10, 18,19, 20, 21};
	bool isLarge = compare_without_anything(d, -1);

	for (int i = 1; i <= 10; ++i)
		cout << fibonacci_n() << endl;


	ARRAY::Array<int> array;
	cout << array;

	int e[] = { 3, 4, 7, 1, 5, 2, 8, 12, 15, 17 };
	const int sizeOfSelected = 10;
	int *f = new int[sizeOfSelected];
	randomSelect(e, f, sizeOfSelected, 10);
	print(f, sizeOfSelected);
	delete[] c;
	delete[] f;

	// Graphic
	/*
	Vertex g[n] = {
		{ 'a', 0, -1, -1},
		{ 'b', 1, -1, -1 },
		{ 'c', 2, -1, -1 },
		{ 'd', 3, -1, -1 },
		{ 'e', 4, -1, -1 },
		{ 'f', 5, -1, -1 },
	};
	int gm[n][n] = {
		{ 0, 1, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 1, 0 },
		{ 0, 0, 0, 0, 1, 1 },
		{ 0, 1, 0, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 0 },
		{ 0, 0, 0, 0, 0, 1 },
	};
	BFS(g, gm, 0);
	DFS(g, gm, 0);
	*/
	//Graphic::G1 g;
	//cout << "The BFS of the graphic is:" << endl;
	//g.BFS(1); cout << endl;
	//cout << "The DFS of the graphic is:" << endl;
	//g.DFS(); cout << endl;
	
	// rescurise
	//char amn[] = "abc";
	//Permutation(amn, amn);
	//int amn[] = { 1, 2, 3 };
	//Amn(amn, 0, 2);

	// DFS
	printf("C53:\n");
	int cmn[] = { 1, 2, 3, 4, 5 };
	std::vector<int> vcmn;
	Cmn(cmn, 3, 0, 4, vcmn);
	/*
	printf("the sum of the number 10 is:\n");
	std::vector<int> sumof;
	SumOf(10, 1, sumof);
	*/
	system("PAUSE");
	return 0;
}
