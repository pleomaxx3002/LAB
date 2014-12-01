// Lab2 1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "LongInt.h"
#include <map>
#define ull unsigned long long
#define ll long long
template <class Type> Type ferma_fact(Type n);
template <typename T> T powmod(T base, T exponent, T &modul);
template <class T> T polard_rho(T n, int  it_count = 1000000);
template <class T> T searchforsqrt(T n);
template <class T> T SQUFOF(T n, int k = 1);
template <class T> T SQ(T &n);
template <typename T> T Lenstra(T n, int k = 10000);
template <typename T> T inversemodule(T, T);
template <typename T> T gcd(T a, T b, T &x, T&y);
template <typename T> T DLogNative(T arg, T base, T module);
template <typename T> T Baby_step_Giant_step(T arg, T base, T module);

template <typename T> class Point_E
{
private:
	T x, y;
	T a, b, n;
public:
	Point_E();
	pair<T, T> Get();
	void Set(T ix, T iy, T ia, T ib, T in);
	Point_E<T> operator*(int);
	Point_E<T> operator+ (Point_E<T>&);
};

template <typename T> pair<T, T> Point_E<T>::Get()
{
	pair<T, T>res;
	res.first = x;
	res.second = y;
	return res;
}
template <typename T> Point_E<T>::Point_E()
{
	x = T(0);
	y = T(0);
	a = b = n = 0;
}


template <typename T> void Point_E<T>::Set(T ix, T iy, T ia, T ib, T in)
{
	x = ix;
	y = iy;
	a = ia;
	b = ib;
	n = in;
}

template <typename T> Point_E<T> Point_E<T>::operator*(int p)
{
	Point_E<T> res;
	Point_E<T> P = (*this);
	res.Set(-1, -1, a, b, n);
	while (p > 0)
	{
		if (p & 1 == 1)
			res = (res.x == -1) ? P : res + P;
		P = P + P;
		p >>= 1;
	}
	return res;
}

template <typename T> Point_E<T> Point_E<T>::operator+(Point_E<T> &B)
{
	if (n <= 0 || a != B.a || b != B.b || n != B.n)
	{
		Point_E <T> c;
		c.Set(-1, -1, -1, -1, -1);
		return c;
	}
	Point_E<T> res = B;
	T lamda;
	if (x == B.x)
	{
		if (y == B.y)
			lamda = (inversemodule((2 * y) % (n), (n)) * ((3 * x * x) % (n)+(a))) % n;
		else
			return Point_E<T>();
	}
	else
	{
		lamda = (B.y - y);
		while (lamda < 0) lamda = lamda + n;
		lamda = (lamda*inversemodule(B.x - x, n)) % n;
	}
	res.x = ((lamda*lamda) % (n) - x);
	while (res.x < 0) res.x = res.x + (n);
	res.x = res.x - B.x;
	while (res.x < 0) res.x = res.x + (n);
	res.y = x - res.x;
	while (res.y < 0) res.y = res.y + (n);
	res.y = (lamda*res.y) % (n) - y;
	while (res.y < 0)res.y = res.y + (n);
	return res;
}


template <typename T> bool RabinMillerPrimeTest(T testednumber, int r = 10);

template <typename T> int JakobiSymbol(T a, T b);

template <typename T> bool Lukas_Selfridge(T number);

template <typename T> bool BPSW(T num);

template <typename T> T phi(T num);

template <typename T> T Generator(T num);

template<typename T> T Polard_Rho_DLog(T base, T num, T module);

template <typename T> vector<pair<T, int>>* Factor(T num);

template <typename T> T Pohlig_Hellman(T a, T b, T num); 

template <typename T> T Index(T a, T b, T num);





int main()
{
	LongInt::LongIntMultMode = Schonhage_Strassen;
	ifstream fIn("input.txt");
	ofstream fOut("output.txt");
	ll m, a, b;
	fIn >> a >> b >> m;
	cout << Pohlig_Hellman(a,b,m) << endl;
	cout << Baby_step_Giant_step(a, b, m) << endl;
	cout << Index(a, b, m) << endl;
	fIn.close();
	fOut.close();
	system("pause");
	return 0;
}










template <class T> T powmod(T b, T e, T &m)
{
	T res = 1;
	while (e > 0)
	{
		if (e % 2 == 1)
			res = (res*b) % m;
		b = (b*b) % m;
		while (b < 0) b = b + m;
		while (res < 0) res = res + m;
		e = e / 2;
	}
	return res;
}

template <class T> T gcd(T a, T b)
{
	if (b < 0) b = -b;
	if (a < 0) a = -a;
	if (a < b)
		swap(a, b);
	while (b != 0)
	{
		a = a%b;
		if (a < b)
			swap(a, b);
	}
	return a;
}

template <class Type> Type ferma_fact(Type n)
{
	Type x = sqrt(n), y = 0, r = x*x - n;
	while (1)
	{
		if (r == 0)
			return (x != y ? (x - y) : (x + y));
		else
		if (r > 0)
		{
			r -= 2 * y + 1;
			++y;
		}
		else
		{
			r += 2 * x + 1;
			++x;
		}
	}
}

template <class T> T polard_rho(T n, int  it_count)
{
	T b0 = rand() % n;
	T b1 = b0;
	T g;
	b1 = (b1*b1) % n;
	if (++b1 == n)
		b1 = 0;
	g = gcd(abs(b1 - b0), n);
	if (g != 1) return g;
	for (int count = 0; count < it_count; count++)
	{
		b0 = (b0*b0) % n;
		if (++b0 == n)
			b0 = 0;
		b1 = (b1*b1) % n;
		++b1;
		b1 = (b1*b1) % n;
		if (++b1 < n)
			b1 = 0;
		g = gcd(abs(b0 - b1), n);
		if (g != 1) return 0;
	}
	return g;
}

template <class T> T searchforsqrt(T n)
{
	T b = sqrt(n);
	for (T i = 2; i < b;++i)
	if (n%i == 0) return i;
	return 1;
}

template <class T> T SQUFOF(T n, int k)
{
	T sq = sqrt(n);
	if (n % 2 == 0) return 2;
	T temp = n - sq*sq;
	if (temp==0) return sq;
	temp = k*n;
	sq = sqrt(temp);
	temp = temp - sq*sq;
	if (temp == 0) return ((temp = gcd(sq, n)) == n) ? 1 : temp;
	T Q1 = 1, P1 = sq, Q2 = temp, r2 = (sq + P1) / Q2;
	T Q0, r1, P0;
	do
	{
		P0 = P1;
		r1 = r2;
		Q0 = Q1; Q1 = Q2;
		P1 = r1*Q1 - P0;
		Q2 = Q0 + (P0 - P1)*r1;
		r2 = (P1 + sq) / Q2;
		temp = sqrt(Q1);
	} while (temp*temp != Q1);
	r1 = (sq + P0) / temp;
	P1 = r1*temp + P0;
	Q1 = temp;
	Q2 = (k*n - P1*P1) / Q1;
	r2 = (sq + P1) / Q2;
	do
	{
		P0 = P1;
		r1 = r2;
		Q0 = Q1, Q1 = Q2;
		P1 = r1*Q1 - P0;
		Q2 = Q0 + (P0 - P1)*r1;
		r2 = (P1 + sq) / Q2;
	} while (P1 != P0);
	temp = gcd(n, P1);
	return (temp==n)? 1 : temp;
}

template <class T> T SQ(T &n)
{
	int i = 1;
	T temp=1;
	while (i <= 1000 && temp == 1)
		temp = SQUFOF(n, i++);
	return temp;
}

template <typename T> T Lenstra(T n, int k)
{
	T x, y, a, b, g, P0, Pxy;
	if (n % 2 == 0)
		return 2;
	if (k > 10000)
		k = 9999;
	ifstream inF("prime.txt");
	vector <int> prime(k);
	for (int i = 0; i < k; i++)
		inF >> prime[i];
	inF.close();
	Point_E<T> P;
	//Initialize
	{
		x = n % ((n - 1) / 3);
		y = n / 3;
		a = 2 * y;
		b = (powmod(y, T(2), n) - ((x*x) % n) - ((a*x) % n)) % n;
		g = gcd(n, 4 * a*a*a + 27 * b*b);
		while (g == n)
		{
			T temp = n - 1;
			x = x + ((rand() % 3) - 1)*(rand() % 200);
			y = y + ((rand() % 3) - 1)*(rand() % 200);
			a = a + ((rand() % 3) - 1)*(rand() % 200);
			if (x > temp) x = temp;
			if (x < 0) x = 0;
			if (y > temp) y = temp;
			if (y<0) y = 0;
			if (a>temp) a = temp;
			if (a<0) a = 0;
			b = (powmod(y, T(2), n) - ((x*x) % n) - ((a*x) % n)) % n;
			g = gcd(n, 4 * a*a*a + 27 * b*b);
		}
		if (g > 1)
			return g;
		P.Set(x, y, a, b, n);
	}
	//Computing
	{
		for (int i = 1; i < k; i++)
		{
			T temp = T(prime[i]), K=T(prime[k-1]);
			while (temp <= K)
			{
				try
				{
					P = P*prime[i];
					temp = temp*prime[i];
				}
				catch (T A)
				{
					return gcd(A, n);
				}
			}
		}
	}

	Point_E<T>P2 = P;
	try
	{
		P2 = P * 2;
		P = P*prime[k - 1];
		while (1)
			P = P + P2;
	}
	catch (T A)
	{
		return gcd(A, n);
	}
}

template <typename T> T gcd(T a, T b, T &x, T &y)
{
	if (a == 0)
	{
		x = 0;
		y = 1;
		return b;
	}
	T x1, y1;
	T d = gcd(b%a, a, x1, y1);
	x = y1 - (b / a)*x1;
	y = x1;
	return d;
}

template <typename T> T inversemodule(T b, T module)
{
	T x, y;
	while (b < 0)
		b = b + module;
	T g = gcd(b, module, x, y);
	if (g != 1)
	{
		throw b;
		return T(0);
	}
	else
	{
		x = x%module;
		if (x < 0)
			x = x + module;
	}
	return x;
}

template <typename T> T DLogNative(T ar, T base, T mod)
{
	T temp = T(1);
	for (T i = T(0); i < mod; ++i)
	{
		if (temp == ar)
			return i;
		temp = (temp*base)%mod;
	}
	return T(-1);
}

template <typename T> T Baby_step_Giant_step(T b, T a, T mod)
{
	if (a == -1)
		return inversemodule(b, mod);
	T m = sqrt(mod) + 1;
	T y = powmod(b, m, mod), am = y;
	map<T, T> mymap;
	for (T i = T(1); i <= m; ++i)
	{
		mymap[y] = i;
		y = (y*am) % mod;
	}
	y = T(1);
	for (T j = 0; j < m; ++j)
	{
		map<T, T>::iterator it;
		it = mymap.find((y*a) % mod);
		if (it != mymap.end())
			return (m*it->second - j);
		y = (y*b) % mod;
	}
}

template <typename T> bool RabinMillerPrimeTest(T m, int r)
{
	T m1 = m - 1, s = 0;
	while (m1>0 && m1 % 2 == 0)
		m1 = m1 / 2, s++;
	srand(time(0));
	for (int i = 0; i < r; i++)
	{
		T a = T(rand()) % (m - 4) + 2;
		for (int j = 0; j < rand() % 10; j++)
			a = (a + T(rand()) % (m - 4)) % (m - 4) + 2;
		T x = powmod(a, m1, m);
		if (x == 1 || x == (m - 1))
			continue;
		for (int j = 0; j < s; j++)
		{
			x = (x*x) % m;
			if (x == 1) return false;
			if (x == m - 1)goto ACycle;
		}
		return false;
	ACycle:;
	}
	return true;
}

template <typename T> int JakobiSymbol(T a, T b)
{
	if (gcd(a, b) != 1) return 0;
	int r = 1;
	if (a < 0)
	{
		a = -a;
		if (b % 4 == 3)
			r = -r;
	}
	while (a!=0)
	{
		int t = 0;
		while (a % 2 == 0)
		{
			a = a / 2;
			++t;
		}
		int temp;
		if ((t & 1) == 1)
		{
			temp = b % 8;
			if (temp == 3 || temp == 5)
				r = -r;
		}
		temp = a % 4;
		if (temp == 3 && temp == b % 4)
			r = -r;
		T c = a;
		a = b%c;
		b = c;
	}
	return r;
}

template <typename T> bool Lukas_Selfridge(T num)
{
	if (num == 2)
		return true;
	if (num % 2 == 0)
		return false;
	T temp = sqrt(num);
	if (temp*temp == num)
		return false;

	//Selfidge
	int d;
	T g;
	for (int tsgn = 1, tabs = 5;; tsgn = -tsgn, tabs = tabs + 2)
	{
		d = tsgn*tabs;
		g = gcd(num, T(tabs));
		if (g > 1)
			return false;
		if (JakobiSymbol(T(d), num) == -1)
			break;
	}
	int p = 1, q = (p*p - d) / 4;

	//Initialization
	g = num + 1;
	T s = 0;
	while (g % 2 == 0)
		g = g / 2, ++s;
	T   u = 0,
		v = 2,
		u2m = 1,
		v2m = p,
		qm = q,
		qm2 = q * 2,
		qkd = q;
	temp = g;

	T inv2 = inversemodule(T(2), num);
	// Lucas
	T u1, v1;
	while (temp > 0)
	{
		if (temp % 2 == 1)
		{
			u1 = u, v1 = v;
			u = ((((u1*v2m) % num + (u2m*v) % num) % num)*inv2) % num;
			v = ((((v1*v2m) % num + (d*u1*u2m) % num) % num)*inv2) % num;
		}
		u1 = u2m, v1 = v2m;
		while (qm2 >= num) qm2 = qm2 - num;
		u2m = (u1*v1) % num;
		v2m = (v1*v1) % num;
		while (v2m < qm2) v2m = v2m + num;
		v2m = v2m - qm2;
		qm = (qm*qm) % num;
		qm2 = qm * 2;
		temp = temp / 2;
	}
	

	if (u == 0 || v == 0)
		return true;

	qm = powmod(T(q), g, num);
	qm2 = 2 * qm; while (qm2 >= num) qm2 = qm2 - num;
	for (int i = 0; i < s; i++)
	{
		v = (v*v) % num;
		while (v < qm2) v = v + num;
		v = v - qm2;
		if (v == 0) return true;
		qm = (qm*qm) % num;
		qm2 = qm * 2; while (qm2 >= num)qm2 = qm2 - num;
	}
	return false;
}

template <typename T> bool BPSW(T num)
{
	int v[1000];
	ifstream inf("prime.txt");
	for (int i = 0; i < 1000; i++)
		inf >> v[i];
	for (int i = 0; i < 1000; i++)
	{
		if (num == v[i])
			return true;
		if (num%v[i] == 0)
			return false;
	}
	
	if (!RabinMillerPrimeTest(num, 100))
		return false;
	return Lukas_Selfridge(num);
}

template <typename T> T phi(T num)
{
	if (num == 1) return 1;
	if (BPSW(num))
		return num - 1;
	T k = SQ(num);
	num = num / k;
	T temp = gcd(num, k);
	return phi(num)*phi(k)*temp / phi(temp);
}

template <typename T> T Generator(T num)
{
	// FACTORIZATION
	T n1 = phi(num), n = n1;
	if (num == 0)
		return 0;
	vector<T> c;
	while (n!=1 && !BPSW(n))
	{
		T temp = SQ(n);
		if (temp == 1) break;
		while (!BPSW(temp))
		{
			temp = temp / SQ(temp);
		}
		while (n%temp == 0)
			n = n / temp;
		c.push_back(temp);
	}
	if (n != 1)
		c.push_back(n);
	int k = c.size();
	srand(time(0));
	T b = 1;

	// Calculating g^(phi(num)/c[i]);
	T g;
	while (b == 1 || b == 0)
	{
		g = rand() % num;
		for (int i = 0; i < rand() % 10; i++)
			g = (g*rand() + rand()) % (num - 1) + 1;
		while (g < 0)g = g + num;
		for (int i = 0; i < k; i++)
		{
			b = powmod(g, n1 / c[i], num);
			if (b == T(1) || b == 0)
				break;
		}
	}
	return g;
}

template <typename T> T Polard_Rho_DLog(T g, T t, T p)
{
	T a = 0, b = 0, x = 1;
	T p1 = phi(p);
	T A = 0, B = 0, X = 1;
	do
	{
		switch (x % 3)
		{
		case T(0): x = (x*x) % p; a = (a * 2) % p1; b = (b * 2) % p1; break;
		case T(1): x = (x*g) % p; a = (a + 1) % p1; break;
		case T(2): x = (x*t) % p; b = (b + 1) % p1; break;
		}

		switch (X % 3)
		{
		case T(0): X = (X*X) % p; A = (A * 2) % p1; B = (B * 2) % p1; break;
		case T(1): X = (X*g) % p; A = (A + 1) % p1; break;
		case T(2): X = (X*t) % p; B = (B + 1) % p1; break;
		}
		switch (X % 3)
		{
		case T(0): X = (X*X) % p; A = (A * 2) % p1; B = (B * 2) % p1; break;
		case T(1): X = (X*g) % p; A = (A + 1) % p1; break;
		case T(2): X = (X*t) % p; B = (B + 1) % p1; break;
		}
	} while (x != X);
	
	A = a - A; while (A < 0) A = A + p1;
	if (A == 0) return -1;
	B = B - b; while (B < 0) B = B + p1;
	T gx, gy;
	X = gcd(B, p1, gx, gy);
	if (X == 1)
		return gx*B;
	b = p1 / X;
	x = (A*gx / X) % p1; while (x < 0) x = x + p1;
	A = powmod(g, x, p);
	a = powmod(g, b, p);
	for (T i = 0; i < X; i++)
	{
		if (A == t)
			return x;

		x = (x + b) % p1;
		A = (A*a) % p;
	}
	return -1;
}

template <typename T> vector<pair<T, int>>* Factor(T num)
{
	vector<pair<T, int>> *f = new vector<pair<T, int>>(0);
	ifstream inf("prime.txt");
	const int N = 1000;
	T temp;
	if (BPSW(num))
	{
		f->push_back(make_pair(num, 1));
		return f;
	}
	for (int i = 0; i < N && num>1; i++)
	{
		inf >> temp;
		int s = 0;
		while (num%temp == 0)
			num = num / temp, ++s;
		if (s != 0)
			f->push_back(make_pair(temp, s));
	}
	if (num == 1)
		return f;
	while (!BPSW(num))
	{
		temp = SQ(num);
		vector<pair<T, int>> *tf = Factor(temp);
		int k = tf->size();
		num = num / temp;
		for (int i = 0; i < k; i++)
		{
			int s = (*tf)[i].second;
			temp = (*tf)[i].first;
			while (num%temp == 0)
				num = num / temp, ++s;
			f->push_back(make_pair(temp, s));
		}
	}
	if (num != 1)
		f->push_back(make_pair(num, 1));
	return f;
}

template <typename T> T Pohlig_Hellman(T a, T b, T num)
{
	T p = phi(num);
	vector<pair<T, int>> *f = Factor(p);
	vector<T>  x(f->size());
	for (int i = 0, n = f->size(); i < n; ++i)
	{
		T q = (*f)[i].first;
		T c = 1, y_1 = 0;
		T a_ = powmod(a, p / q, num), b_;
		vector<T> y((*f)[i].second);
		for (T k = 0; k < (*f)[i].second; ++k)
		{

			if (k == 0)
			{
				b_ = powmod(b, p / q, num);
				y[0] = (b_ == 1) ? 0 : ((b_ == a_) ? 1 : Baby_step_Giant_step(a_, b_, num));
			}
			else
			{
				T temp = powmod(q, k - 1, num);
				c = (c*powmod(a, y[int(k) - 1] * temp, num)) % num;
				b_ = powmod(b*inversemodule(c, num), p / ((temp*q*q) % num), num);
				y[k] = (b_ == 1) ? 0 : ((b_ == a_) ? 1 : Baby_step_Giant_step(a_, b_, num));// mod=num. mod=powmod(q,k+1,num);
			}
		}
		b_ = 1;
		x[i] = 0;
		(*f)[i].first = powmod((*f)[i].first, (T)(*f)[i].second, num);
		for (int j = 0; j < (*f)[i].second; j++)
			x[i] = (x[i] + y[j] * b_) % (*f)[i].first, b_ = b_*q;
	}
	vector<T>xi(x);
	for (int i = 0, n = f->size(); i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			xi[i] = inversemodule((*f)[j].first, (*f)[i].first)*(xi[i] - xi[j]);
			while (xi[i] < 0) xi[i] = xi[i] + (*f)[i].first;
			xi[i] = xi[i] % (*f)[i].first;
		}
	}
	T res = 0;
	for (int i = 0, n = f->size(); i < n; i++)
	{
		T temp = xi[i];
		for (int j = 0; j < i; j++)
			temp = temp*(*f)[j].first;
		res = res + temp;
	}
	return res;
}

template <typename T> T Index(T a, T b, T num)
{
	int c = 2;
	T n1 = phi(num);
	ifstream inF("prime.txt");
	vector<int> prime;
	double L = exp(log(num)*log(log(num))) + 1;
	if (L > num) L = num;
	prime.resize(int(L));
	int i;
	for (i = 0;; ++i)
	{
		inF >> prime[i];
		if (prime[i] > L) break;
	}
	prime.resize(i);
	int N = i + 1;
	int n = 0;
	T ak;
	srand(time(0));
	T k = 1;
	vector<vector<T>> m(N + c, vector<T>(N + 1, 0));
	vector<T> m1(N + 1, 0);
	while (true)
	{
		m.resize(N + c, vector<T>(N + 1, 0));;
		//Random k: a^k=П(prime[i]);
		while (n < N + c)
		{
			
			k = (k + powmod(T(rand() % 234), T(rand() % 10), n1)) % (n1 - 1) + 1;
			ak = powmod(a, k, num);
			if (ak == b)return k;
			if (ak == 1) continue;
			for (int i = 0; i < N, ak>1; i++)
			{
				while (ak%prime[i] == 0)
					ak = ak / prime[i], ++m[n][i];
			}
			if (ak != 1) m[n] = vector<T>(N + 1, 0);
			else
				m[n][N] = k, n++;
		}

		//Lineal System
		for (int i = 0; i < N; i++)
		{
			int j = i;
			while (j < n && (m[j][i] == 0 || gcd(m[j][i], n1) != 1))
				j++;
			if (j == n) goto again;
			swap(m[i], m[j]);
			T X, Y;
			for (int I = 0; I < n; I++)
			{ 
				X = inversemodule(m[I][i], n1);
				if (I == i)
					for (int J = i; J <= N; J++)
						m[I][J] = (m[I][J] * X) % (n1);
				else
				{
					X = -((X*m[I][i]) % (n1));
					while (X < 0) X = X + (n1);
					for (int J = i; J <= N; ++J)
					{
						m[I][J] = (m[I][J] + X*m[j][J]) % (n1);
						while (m[I][J] < 0) m[I][J] = m[I][J] + (n1);
					}
				}
			}

		}
		//Random k: b*a^k
		k = 1;
		bool F = true;
		while (F)
		{
			while (k != 0)
			{
				k = (k + powmod(T(rand() % 234), T(rand() % 10), (n1))) % (num - 2) + 1;
				if (k != 0) ak = powmod(a, k, num);
				if (ak == b) return k;
			}
			ak = b*ak;
			if (ak == 1) continue;
			for (int i = 0; i < N, ak>1; i++)
			{
				while (ak%prime[i] == 0)
					ak = ak / prime[i], ++m1[i];
			}
			if (ak != 1) m1 = vector<T>(N + 1, 0);
			else
				F=false;
		}
		T res = (n - 1) - k;;
		for (int i = 0; i < N; i++)
			res = (res + m1[i] * m[i][N]) % (n - 1);
		return res;

	again:;
		c++;
	}
}


