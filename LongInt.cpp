#include "stdafx.h"
#include "LongInt.h"


LongInt (*LongInt::LongIntMultMode)(LongInt &A, LongInt &B)=NULL;
//constructor&destructor
LongInt::LongInt()
{
	sgn = 1;
	dot = -1;
}

LongInt::LongInt(int N)
{
	sgn = (N < 0) ? -1 : 1;
	dot = -1;
	N = abs(N);
	int i = 0;

	while (N>0)
	{
		num.push_back(N%base);
		N /= base;
		i++;
	}
}

LongInt::LongInt(LongInt &A)
{
	num = A.num;
	sgn = A.sgn;
	dot = A.dot;
}

LongInt::LongInt(vector<int>&A)
{
	num = A;
	sgn = 1;
	dot = -1;
}

LongInt::~LongInt()
{
}



LongInt abs(LongInt &B)
{
	LongInt res(B);
	res.sgn = abs(B.sgn);
	return res;
}


//comparisoin
bool LongInt::operator> (LongInt &B)
{
	if (sgn > B.sgn)
		return true;
	if (sgn < B.sgn)
		return false;
	if (sgn < 0)
	{
		Normalize();
		B.Normalize();
		if (num.size()>B.num.size())
			return false;
		else
		if (num.size() < B.num.size())
			return true;
		for (int i = num.size() - 1; i >= 0; i--)
		{
			if (num[i] < B.num[i])
				return true;
			else
			if (num[i] > B.num[i])
				return false;
		}
		return false;
	}
	else
	{
		Normalize();
		B.Normalize();
		if (num.size() < B.num.size())
			return false;
		else
		if (num.size() > B.num.size())
			return true;
		for (int i = num.size() - 1; i >= 0; i--)
		{
			if (num[i] > B.num[i])
				return true;
			else
			if (num[i] < B.num[i])
				return false;
		}
		return false;
	}
}

bool LongInt::operator>=(LongInt &B)
{
	{
		if (sgn > B.sgn)
		return true;
		if (sgn < B.sgn)
			return false;
		if (sgn < 0)
		{
			Normalize();
			B.Normalize();
			if (num.size()>B.num.size())
				return false;
			else
			if (num.size() < B.num.size())
				return true;
			for (int i = num.size() - 1; i >= 0; i--)
			{
				if (num[i] < B.num[i])
					return true;
				else
				if (num[i] > B.num[i])
					return false;
			}
			return true;
		}
		else
		{
			Normalize();
			B.Normalize();
			if (num.size() < B.num.size())
				return false;
			else
			if (num.size() > B.num.size())
				return true;
			for (int i = num.size() - 1; i >= 0; i--)
			{
				if (num[i] > B.num[i])
					return true;
				else
				if (num[i] < B.num[i])
					return false;
			}
			return true;
		}
	}
}

bool LongInt::operator== (LongInt &B)
{
	if (sgn != B.sgn)
		return false;
	if (num.size() != B.num.size())
		return false;
	for (int i = 0; i < num.size();i++)
	if (B.num[i] != num[i])
		return false;
	return true;
}

bool LongInt::operator<(LongInt &B)
{
	return !((*this) >= B);
}

bool LongInt::operator<=(LongInt &B)
{
	return !((*this) >B);
}

bool LongInt::operator!=(LongInt &B)
{
	return !((*this) == B);
}


//+-
LongInt  LongInt::operator+(int a)
{
	return (*this) + LongInt(a);
}

LongInt& LongInt::operator++()
{
	(*this) = (*this) + 1;
	return (*this);
}

LongInt operator++(LongInt &A, int B)
{
	LongInt old(A);
	++A;
	return old;
}

LongInt  LongInt::operator+ (LongInt &B)
{
	if (sgn != B.sgn)
	{
		if (abs((*this)) > abs(B))
		{
			if (sgn > 0)
				return (*this) - B.Minus();
			return (this->Minus() - B).Minus();
		}
		else
		{
			if (sgn > 0)
				return (B.Minus() - (*this)).Minus();
			return B - (*this).Minus();
		}
	}
	LongInt res(*this);
	res.num.resize(max(B.num.size(), num.size()) * 2, 0);
	int minimum = B.num.size();
	int add = 0;
	for (int i = 0; i < minimum; i++)
	{
		res.num[i] += B.num[i] + add;
		add = res.num[i] / base;
		res.num[i] %= base;
	}
	while (add != 0)
	{
		res.num.push_back(0);
		res.num[minimum] += add;
		add = res.num[minimum] / base;
		res.num[minimum] %= base;
		minimum++;
	}
	res.Normalize();
	return (res);
}

LongInt  operator +(int a, LongInt &B)
{
	return  B + a;
}

LongInt  LongInt::operator-(LongInt &B)
{
	if (sgn != 1 || B.sgn != 1)
		return (*this) + B.Minus();
	if ((*this) < B)
		return (B - (*this)).Minus();
	int n = B.num.size();
	LongInt res(*this);
	int add = 0;
	for (int i = 0; i < n; i++)
	{
		res.num[i] -= B.num[i]+add;
		add = 0;
		if (res.num[i] < 0)
			res.num[i] += base, add = 1;
	}
	while (add != 0)
	{
		res.num[n] -= add;
		add = 0;
		if (res.num[n] < 0)
			res.num[n] += base, add = 1, n++;
	}
	res.Normalize();
	return res;
}

LongInt LongInt::operator-(int B)
{
	return (*this) - LongInt(B);
}



LongInt LongInt::Minus()
{
	LongInt res(*this);
	res.sgn = -sgn;
	return res;
}

void LongInt::Normalize()
{
	int i; 
	for (i = num.size() - 1; i >= 0 && num[i]==0; i--);
	num.resize(i + 1);
}


LongInt LongInt::operator* (int B)
{
	ull temp = 0;
	LongInt res(*this);
	if (B < 0)
		res.sgn = -res.sgn;
	for (int i = 0; i < num.size(); i++)
	{
		temp = ull(B)*num[i]+temp;
		res.num[i] = temp % base;
		temp /= base;
	}
	if (temp != 0)
		res.num.push_back(temp);
	return res;
}

LongInt operator* (int A, LongInt &B)
{
	return B*A;
}

//I/O

istream & operator >> (istream &stream, LongInt &A)
{
	string h;
	stream >> h;
	if (h[0] == '-')
		A.sgn = -1, reverse(h.begin(),h.end()), h.pop_back();
	else
	if (h[0] == '+')
		A.sgn = 1, reverse(h.begin(), h.end()), h.pop_back();
	else
		A.sgn = 1, reverse(h.begin(), h.end());
	int numsize = log10(LongInt::base);
	int hsize = h.size()/numsize;
	A.num.assign(hsize + 1, 0);
	for (int i = 0; i < hsize; i++)
	{
		int t = 1;
		for (int j = 0; j < numsize; j++)
		{
			if (h[i*numsize + j]<'0' || h[i*numsize + j]>'9')
			{
				A.num.assign(0, 0);
				cout << "Input error\n";
				return stream;
			}
			A.num[i] += (h[i*numsize + j]-'0') * t;
			t *= 10;
		}
		if (A.num[i]>A.base)
		{
			A.num.assign(0,0);
			cout << "Input error\n";
			return stream;
		}
	}
	int t = 1;
	for (int i = hsize*numsize; i < h.size(); i++)
	{
		if (h[i]<'0' || h[i]>'9')
		{
			A.num.assign(0, 0);
			cout << "Input error\n";
			return stream;
		}
		A.num[hsize] += t*(h[i] - '0');
		t *= 10;
	}
	A.Normalize();
	return stream;
}

ostream & operator << (ostream& stream, const LongInt& A)
{
	if (A.num.size() == 0)
	{
		stream << 0;
		return stream;
	}
	if (A.sgn == -1)
		stream << '-';
	int i = A.num.size() - 1;
	stream << A.num[i--];
	int numsize = log10(A.base);
	for (; i >= 0; i--)
		stream.fill('0'), stream.width(numsize), stream << A.num[i];
	stream << endl;
	return stream;
}





LongInt Karatsuba(LongInt &A, LongInt &B)
{
	int k;
	for (k = A.num.size() - 1; k >= 0; k--)
	if (A.num[k] != 0)
		break;
	A.num.resize(k + 1);
	for (k = B.num.size() - 1; k >= 0; k--)
	if (B.num[k] != 0)
		break;
	B.num.resize(k + 1);
	int s1 = A.num.size();
	int s2 = B.num.size();
	LongInt res;
	if (s1 == 0)
		return LongInt();
	else
	if (s2 == 0)
		return LongInt();
	if (s1 == 1)
		return B*A.num[0];
	else
	if (s2 == 1)
		return A*B.num[0];
	if (s1 < s2)
		A.num.resize(s2, 0), s1 = s2;
	if (s2 < s1)
		B.num.resize(s1, 0), s2 = s1;
	s2 = s1 / 2;
	res.num.assign(s1,0);
	LongInt L1(A);
	L1.num.resize(s2);
	L1.sgn = 1;
	LongInt R1 = A >> s2;
	R1.sgn = 1;
	LongInt L2(B);
	L2.sgn = 1;
	L2.num.resize(s2);
	LongInt R2 = B >> s2;
	R2.sgn = 1;
	LongInt z0 = Karatsuba(L1, L2);
	LongInt z1 = Karatsuba((L1 + R1), (L2 + R2));
	LongInt z2 = Karatsuba(R1, R2);
	if (A.sgn*B.sgn > 0)
		return (z2 << (s2*2)) + ((z1 - z2 - z0) << (s2)) + z0;
	else
		return ((z2 << (s2*2)) + ((z1 - z2 - z0) << (s2)) + z0).Minus();
}

LongInt LongInt::operator<< (int a)
{
	LongInt res = LongInt(*this);
	res.num.resize(num.size() + a);
	for (int i = res.num.size() - 1; i >= a; i--)
		res.num[i] = res.num[i - a];
	for (int i = 0; i < a; i++)
		res.num[i] = 0;
	return res;
}

LongInt LongInt::operator>> (int a)
{
	LongInt res = LongInt(*this);
	if (a >= res.num.size())
		return LongInt(0);
	for (int i = a; i < res.num.size(); i++)
	{
		res.num[i - a] = res.num[i];
	}
	res.num.resize(num.size() - a);
	return res;
}

LongInt LongInt::operator*(LongInt &B)
{
	clock_t start=clock();
	LongInt res(LongIntMultMode(*this, B));
	cout << "Multiplication duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
	return res;
}

LongInt Toom(LongInt &A, LongInt &B)
{
	A.Normalize(), B.Normalize();
	int s1 = A.num.size();
	int s2 = B.num.size();
	LongInt res;
	if (s1 < 3 || s2 < 3)
		return Native(A, B);
	s1 = max(s1, s2);
	A.num.resize(s1);
	B.num.resize(s1);
	s2 = (s1 % 3 == 0) ? (s1 / 3) : s1 / 3 + 1;
	LongInt M0=abs(A);
	LongInt M1(M0), M2(M0);
	M0.num.resize(s2);
	M1.num.erase(M1.num.begin(), M1.num.begin() + s2); M1.num.erase(M1.num.begin() + s2, M1.num.end());
	M2.num.erase(M2.num.begin(), M2.num.begin() + 2 * s2);
	LongInt N0(abs(B)), N1(N0), N2(N0);
	N0.num.resize(s2);
	N1.num.erase(N1.num.begin(), N1.num.begin() + s2); N1.num.erase(N1.num.begin() + s2, N1.num.end());
	N2.num.erase(N2.num.begin(), N2.num.begin() + 2 * s2);
	LongInt R0 = Toom(M0, N0);
	LongInt R1 = Toom(M0 + M1 + M2, N0 + N1 + N2);
	LongInt R_1 = Toom(N0 - N1 + N2, M0 - M1 + M2);
	LongInt R_2 = Toom(4 * N2 - 2 * N1 + N0, 4 * M2 - 2 * M1 + M0);
	LongInt RInf = Toom(M2, N2);
	LongInt r3 = (R_2 - R1) / 3;
	LongInt r1 = (R1 - R_1) / 2;
	LongInt r2 = R_1 - R0;
	r3 = (r2 - r3) / 2 + 2*RInf;
	r2 = r2 + r1 - RInf;
	r1 = r1 - r3;
	res = R0 + (r1 << s2) + (r2 << (2 * s2)) + (r3 << (3 * s2)) + (RInf << (4 * s2));
	if (A.sgn*B.sgn > 0)
		return res;
	else
		return res.Minus();
}

LongInt LongInt::operator/(int B)
{
	LongInt res(*this);
	int b;
	b = abs(B);
	ll carry = 0;
	for (int i = num.size() - 1; i >= 0; i--)
	{
		ll cur = res.num[i] + carry*base;
		res.num[i] = cur / b;
		carry = cur%b;
	}
	res.Normalize();
	if (B < 0) res.sgn = -res.sgn;
	return res;
}

LongInt Native(LongInt&A, LongInt&B)
{
	LongInt res,temp;
	res.num.resize(max(A.num.size(), B.num.size()) * 2, 0);
	temp.num.resize(res.num.size(), 0);
	for (int i = 0; i < A.num.size(); i++)
	{
		ll ad = 0, cur=0;
		for (int j = 0; j < B.num.size(); j++)
		{
			ad = ll(A.num[i])*B.num[j] + cur;
			temp.num[j] = ad % A.base;
			cur = ad/A.base;
		}
		temp.num[B.num.size()] = cur;
		res = res + (temp << i);
		temp.num.assign(temp.num.size(), 0);
	}
	res.sgn = A.sgn*B.sgn;
	return res;
}

void LongInt::operator=(const LongInt& B)
{
	if (this == &B) return;
	num = B.num;
	sgn = B.sgn;
	return;
}

void FFT(vector< complex <double> > &a, bool dir)
{
	int n = a.size();
	if (n <= 1)
		return;
	vector< complex <double> > a0(n / 2), a1(n / 2);
	for (int i = 0, j = 0; i < n; i += 2, j++)
	{
		a0[j] = a[i];
		a1[j] = a[i + 1];
	}
	FFT(a0, dir);
	FFT(a1, dir);
	double ang = 2 * M_PI / n * (dir ? 1 : -1);
	complex<double> w(1.0), w_n(cos(ang), sin(ang));
	n /= 2;
	for (int i = 0; i < n; i++)
	{
		a[i] = a0[i] + w*a1[i];
		a[i + n] = a0[i] - w*a1[i];
		if (!dir)
			a[i] /= 2, a[i + n] /= 2;
		w *= w_n;
	}
}

LongInt LongInt::InverseCook()
{
	Normalize();
	if (num.size() == 0)
		return 0;
	double t = LongInt::base / double(num[num.size()-1]);
	LongInt z;
	if (num[num.size() - 1] == 1)
		z.num.push_back(0), z.num.push_back(1);
	else
		z.num.push_back(int(t*base) % base), z.num.push_back(floor(t));
	z.dot = num.size() + 1;
	int k = 2;
	int n = 3; 
	LongInt temp;

	for (int i = 0; i < 10 || z.dot < 20*num.size(); i++)
	{
		temp = Schonhage_Strassen(Schonhage_Strassen(*this, z), z);
		z = ((2 * z) << (z.dot)) - (temp);
		z.dot = z.dot*2;
		//z.dot -= z.num.size() - num.size() - 1;
		//z = z >> (z.num.size() - num.size() - 1);

	}
	return z;
}

LongInt Schonhage_Strassen(LongInt &A, LongInt &B)
{
	vector< complex<double> > a(A.num.begin(), A.num.end()), b(B.num.begin(), B.num.end());
	int n = 1;
	int temp = max(a.size(), b.size());
	while (n < temp) n = n << 1;
	n = n << 1;
	a.resize(n, 0);
	b.resize(n, 0);
	FFT(a, true);
	FFT(b, true);
	for (int i = 0; i < n; i++)
		a[i] *= b[i];
	FFT(a, false);
	LongInt res(A);
	res.num.resize(n);
	ll carry = 0;
	for (int i = 0; i < n; i++)
	{
		carry += a[i].real() + 0.5;
		res.num[i] = carry%LongInt::base;
		carry /= LongInt::base;
	}
	res.Normalize();
	return res;
}

LongInt LongInt::operator/ (LongInt &B)
{
	clock_t start = clock();
	LongInt res;
	if (B.num.size() == 0)
		return 0;
	LongInt c = B.InverseCook();
	res = Schonhage_Strassen((*this), c) >> (c.dot);
	cout << "Dividing duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
	return res;
}

LongInt LongInt::operator% (LongInt &B)
{
	B.Normalize();
	int sg = sgn;
	sgn = abs(sgn);
	if (B.num.size() == 0)
		return 0;
	LongInt c = (*this) - Schonhage_Strassen((*this) / B, B);
	while (c < LongInt(0))
		c = c + B;
	while (c>=B)
		c = c - B;
	sgn = sg;
	return c;
}


LongInt Schonhage(LongInt &A, LongInt &B)
{
	int q = 2, qp = 1;
	int p = 18 * q + 8;
	int s = 2 * max(A.num.size(), B.num.size());
	while (p < s)
	{
		q = 3 * q - 1;
		p = 18 * q + 8;
	}
	if (A.num.size() <= (6 * q + 8) || B.num.size() <= (6 * q + 8))
		return Native(A, B);
	if (p <= 44)
		return Native(A, B);
	LongInt m[6];
	m[0] = (LongInt(LongInt::base) << (6 * q - 1)) - LongInt(LongInt::base - 1);
	m[1] = (LongInt(LongInt::base) << (6 * q + 1)) - LongInt(LongInt::base - 1);
	m[2] = (LongInt(LongInt::base) << (6 * q + 2)) - LongInt(LongInt::base - 1);
	m[3] = (LongInt(LongInt::base) << (6 * q + 3)) - LongInt(LongInt::base - 1);
	m[4] = (LongInt(LongInt::base) << (6 * q + 5)) - LongInt(LongInt::base - 1);
	m[5] = (LongInt(LongInt::base) << (6 * q + 7)) - LongInt(LongInt::base - 1);
	LongInt u[6];
	for (int i = 0; i < 6; i++)
		u[i] = A%m[i];
	LongInt v[6];
	for (int i = 0; i < 6; i++)
		v[i] = B%m[i];
	for (int i = 0; i < 6; i++)
		u[i] = (Schonhage(u[i], v[i]) % m[i]);
	v[0] = u[0];
	for (int i = 1; i < 6; i++)
	{
		LongInt c(u[i]);
		for (int j = 0; j < i; j++)
			c = Schonhage((c - v[j]), m[j].pow(m[i] - LongInt(2), m[i]) % m[i]) % m[i];
		v[i] = c;
	}
	return Schonhage((Schonhage((Schonhage((Schonhage((Schonhage(v[5], m[4]) + v[4]), m[3]) + v[3]), m[2]) + v[2]), m[1]) + v[1]), m[0]) + v[0];
}

LongInt LongInt::pow(LongInt exp, LongInt &mod)
{
	LongInt res(1), a(*this), inv(mod.InverseCook());
	while (exp.num.size()>0)
	{
		if (exp.num[0] % 2 == 1)
		{
			res = Schonhage_Strassen(res, a);
			res = res - Schonhage_Strassen(Schonhage_Strassen(res, inv) >> inv.dot, mod);
			while (res < LongInt(0))
				res = res + mod;
			while (res>= mod)
				res = res - mod;

		}
		a = Schonhage_Strassen(a, a);
		a = a - Schonhage_Strassen((Schonhage_Strassen(a, inv) >> inv.dot),mod);
		while (a < LongInt(0))
			a = a + mod;
		while (a >= mod)
			a = a - mod;
		exp = exp / 2;
		exp.Normalize();
	}
	res.Normalize();
	return res;
}

LongInt gcd(LongInt A, LongInt B)
{
	if (A < B)
	{
		LongInt temp(A);
		A = B;
		B = temp;
	}
	while (B.num.size() > 0)
	{
		A = A%B;
		A.Normalize();
		if (A < B)
		{
			LongInt temp(A);
			A = B;
			B = temp;
		}
	}
	A.Normalize();
	return A;
}

bool LongInt::PrimeLukasLehmar()
{
	clock_t start = clock();
	LongInt L(4);
	LongInt M(LongInt(2).pow((*this)) - 1);
	LongInt InvM = M.InverseCook();
	(*this) = (*this) - LongInt(2);
	for (LongInt i(0); i < (*this); i = i + 1)
	{		
		L = L*L - 2;
		L = L - (((L)*InvM) >> InvM.dot)*M;
		while (L < LongInt(0))
			L = L + M;
		while (L >= M)
			L = L - M;
	}
	(*this) = (*this) + 2;
	if (L.num.size() == 0)
	{
		cout << "M_q is prime" << endl;
		cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
		return true;
	}
	cout << "M_q is composite" << endl;
	cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
	return false;
}

LongInt LongInt::pow(LongInt E)
{
	LongInt res(1);
	LongInt A(*this);
	while (!E.num.empty())
	{
		if (E.num[0] % 2 == 1)
			res = res*A;
		A = A*A;
		E = E / 2;
	}
	return res;
}

bool LongInt::PrimeMillerRabin(LongInt &r)
{
	clock_t start = clock();
	if ((*this) <= LongInt(3))
	{
		if ((*this) <= LongInt(0))
		{
			cout << "Input data is incorrect: input value is not positive" << endl;
			return false;
		}
		if ((*this) == LongInt(1))
		{
			cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
			return false;
		}
		cout << "Prime" << endl;
		cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
		return true;
	}
	if (mod2() == 0)
	{
		cout << "Composite" << endl;
		cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
	}

	(*this) = (*this) - LongInt(1);
	LongInt d(*this);
	LongInt s(0);
	LongInt a;
	LongInt x;
	while (d.num[0] % 2 == 0)
		d = d / 2, ++s;
	srand(time(NULL));
	for (LongInt i = 0; i < r; ++i)
	{
		//genereting random a
		a.num.resize(rand() % (num.size()-1)+1);
		if (num[num.size() - 1] == 1 && a.num.size() == num.size())
			a.num.pop_back();
		for (int i = a.num.size()-2; i >=0; --i)
			a.num[i] = rand() % base;
		a.num[a.num.size() - 1] = rand() % (base - 1) + 1;
		if (a == LongInt(1))
			++a;


		x = a.pow(d, (*this) + 1);
		if (x == LongInt(1) || x == (*this)) continue;
		for (LongInt i = 1; i < s; ++i)
		{
			x = (x*x) % ((*this) + 1);
			if (x == LongInt(1))
			{
				++(*this);
				cout << "Composite" << endl;
				cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
			}
			if (x == (*this))
				goto WitnessLoop;
		}
		cout << "Composite" << endl;
		++(*this);
		cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
		return false;
	WitnessLoop:;
	}
	cout << "Probably prime" << endl;
	++(*this);
	cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
	return true;
}

bool LongInt::PrimeSolovayStrassen(LongInt &r)
{
	clock_t start = clock();
	if ((*this) <= LongInt(3))
	{
		if ((*this) <= LongInt(0))
		{
			cout << "Input data is incorrect: input value is not positive" << endl;
			return false;
		}
		if ((*this) == LongInt(1))
		{
			cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
			return false;
		}
		cout << "Prime" << endl;
		cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
		return true;
	}
	if (mod2() == 0)
	{
		cout << "Composite" << endl;
		cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
		return false;
	}
	LongInt a;
	(*this) = (*this) - 1;
	for (LongInt k = 0; k<r; k++)
	{
		//random a generating
		a.num.resize(rand() % (num.size() - 1) + 1);
		if (num[num.size() - 1] == 1 && a.num.size() == num.size())
			a.num.pop_back();
		for (int i = a.num.size() - 2; i >= 0; --i)
			a.num[i] = rand() % base;
		a.num[a.num.size() - 1] = rand() % (base - 1) + 1;
		if (a == LongInt(1))
			++a;
		if (gcd(a, (*this) + 1) > LongInt(1))
		{
			cout << "Composite" << endl;
			++(*this);
			cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
			return false;
		}
		LongInt temp=(a.pow((*this) / 2, (*this) + 1) - JakobiSymbol(a, (*this) + 1));
		temp.Normalize();
		if (temp != ((*this)+1) && !temp.num.empty())
		{
			cout << "Composite" << endl;
			++(*this);
			cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
			return false;
		}

	}
	cout << "Probably Prime" << endl;
	++(*this);
	cout << "Counting duration:" << (clock() - start) / double(CLOCKS_PER_SEC) << "sec" << endl;
	return true;
}


int JakobiSymbol(LongInt A, LongInt B)
{
	if (gcd(abs(A), abs(B)) > LongInt(1))
		return 0;
	int res = 1;
	if (B.num.empty())
	{
		if (A.num.size() == 1 && A.num[0] == 1)
			return 1;
		else
			return 0;
	}
	if (A.sgn < 0)
	{
		A.sgn = -A.sgn;
		if (LongInt::base>100)
		{
			if (B.num[0] % 4 == 3)
				res = -res;
		}
		else
		{
			int t = B.num[0];
			if (B.num.size() > 1)
				t +=LongInt::base*B.num[1];
			if (t % 4 == 3)
				res = -res;
		}
	}
	do
	{
		LongInt t = 0;
		while (!A.num.empty() && A.num[0] % 2 == 0)
		{
			++t;
			A = A / 2;
		}
		if (t.mod2() == 1)
		{
			LongInt temp = B%LongInt(8);
			if (temp.num[0] == 3 || temp.num[0] == 5)
				res = -res;
		}
		if (A.mod4() == 3 && B.mod4() == 3)
			res = -res;
		LongInt c = A;
		A = B%c;
		B = c;
	} while (A.num.size() != 0);
	return res;
}

inline int LongInt::mod4()
{
	if (num.size() == 0)
		return 0;
	if (num.size() == 1)
		return num[0] % 4;
	return (num[0] + ll(base)*ll(num[1])) % 4;
}

inline int LongInt::mod2()
{
	if (num.size() == 0)
		return 0;
	return num[0] % 2;
}

LongInt sqrt(LongInt &A)
{
	LongInt mi(LongInt(1)<<(A.num.size()/2)), ma(LongInt(1)<<(A.num.size()/2 + (A.num.size()%2 ? 2 : 1))), temp;
	while (mi < ma)
	{
		temp = ((mi + ma) / 2).pow(2);
		if (temp == A)
			return temp;
		if (temp < A)
			mi = (mi + ma) / 2;
		else
			ma = (mi + ma) / 2;
	}
	if (ma.pow(2) == A)
		return ma;
	ma = ma << 1;
	ma.num[0] = 1;
	ma.dot++;
	return ma;
}