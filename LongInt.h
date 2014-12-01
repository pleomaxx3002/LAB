#pragma once
class LongInt
{
public:
	LongInt();
	LongInt(int);
	LongInt(LongInt&);
	LongInt(vector<int>&);
	~LongInt();
	LongInt  operator + (LongInt &);
	LongInt  operator + (int);
	LongInt operator-(int);
	friend LongInt operator++(LongInt&, int);
	LongInt &operator++();
	friend LongInt (operator +) (int, LongInt &);
	LongInt operator- (LongInt&);
	bool operator > (LongInt &);
	bool operator == (LongInt &);
	bool operator >= (LongInt &);
	bool operator < (LongInt &);
	bool operator <= (LongInt &);
	bool operator!= (LongInt &);
	friend LongInt  abs(LongInt &);
	friend LongInt sqrt(LongInt&);
	LongInt Minus();


	friend istream & operator >>(istream&, LongInt&);
	friend ostream & operator << (ostream&, const LongInt&);


	void operator=(const LongInt&);
	LongInt operator*(LongInt&);
	LongInt operator * (int);
	LongInt operator << (int);
	LongInt operator >> (int);
	LongInt operator/(int);
	LongInt operator/ (LongInt &);
	LongInt operator% (LongInt &);

	LongInt pow(LongInt, LongInt&);
	static LongInt(*LongIntMultMode)(LongInt &A, LongInt &B);
	friend LongInt operator*(int, LongInt&);
	friend LongInt Schonhage(LongInt &, LongInt &);
	friend LongInt Schonhage_Strassen(LongInt&, LongInt&);
	friend LongInt Native(LongInt&, LongInt&);
	friend LongInt Karatsuba(LongInt&, LongInt&);
	friend LongInt Toom(LongInt&, LongInt &);
	LongInt InverseCook();
	friend LongInt gcd(LongInt, LongInt);
	LongInt pow(LongInt);
	bool PrimeLukasLehmar();
	bool PrimeMillerRabin(LongInt &);
	bool PrimeSolovayStrassen(LongInt &);
	friend int JakobiSymbol(LongInt, LongInt);
	inline int mod4();
	inline int mod2();
private:
	static const int base = 10;
	vector <int> num;
	void Normalize();
	int sgn, dot;
};