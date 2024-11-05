/*
complex.h
*/

#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#if defined(__GNUC__) && !defined(__NEC__) && !defined(__FUJITSU) && !defined(__CLANG_FUJITSU)
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <math.h>

static inline d_complex_t d_complex(double r, double i)
{
	d_complex_t z;

	z.r = r;
	z.i = i;

	return z;
}

static inline d_complex_t d_add(d_complex_t a, d_complex_t b)
{
	return d_complex(a.r + b.r, a.i + b.i);
}

static inline d_complex_t d_sub(d_complex_t a, d_complex_t b)
{
	return d_complex(a.r - b.r, a.i - b.i);
}

static inline d_complex_t d_mul(d_complex_t a, d_complex_t b)
{
	return d_complex((a.r * b.r) - (a.i * b.i), (a.r * b.i) + (a.i * b.r));
}

static inline d_complex_t d_div(d_complex_t a, d_complex_t b)
{
	if ((fabs(b.r) <= 0) && (fabs(b.i) <= 0)) return d_complex(0, 0);
	return d_complex(((a.r * b.r) + (a.i * b.i)) / ((b.r * b.r) + (b.i * b.i)),
	                 ((a.i * b.r) - (a.r * b.i)) / ((b.r * b.r) + (b.i * b.i)));
}

static inline d_complex_t d_rmul(double r, d_complex_t z)
{
	return d_complex(r * z.r, r * z.i);
}

static inline d_complex_t d_add3(d_complex_t a, d_complex_t b, d_complex_t c)
{
	return d_complex(a.r + b.r + c.r, a.i + b.i + c.i);
}

static inline d_complex_t d_add4(d_complex_t a, d_complex_t b, d_complex_t c, d_complex_t d)
{
	return d_complex(a.r + b.r + c.r + d.r, a.i + b.i + c.i + d.i);
}

static inline double d_abs(d_complex_t z)
{
	return sqrt((z.r * z.r) + (z.i * z.i));
}

static inline double d_norm(d_complex_t z)
{
	return ((z.r * z.r) + (z.i * z.i));
}

static inline double d_rad(d_complex_t z)
{
	return atan2(z.i, z.r);
}

static inline double d_deg(d_complex_t z)
{
	return atan2(z.i, z.r) * 45.0 / atan(1.0);
}

static inline d_complex_t d_exp(double x)
{
	return d_complex(cos(x), sin(x));
}

static inline d_complex_t d_cos(d_complex_t z)
{
	return d_add(d_rmul(0.5 * exp(-z.i), d_complex(+cos(z.r), +sin(z.r))),
	             d_rmul(0.5 * exp(+z.i), d_complex(+cos(z.r), -sin(z.r))));
}

static inline d_complex_t d_sin(d_complex_t z)
{
	return d_sub(d_rmul(0.5 * exp(-z.i), d_complex(+sin(z.r), -cos(z.r))),
	             d_rmul(0.5 * exp(+z.i), d_complex(-sin(z.r), -cos(z.r))));
}

static inline d_complex_t d_sqrt(d_complex_t z)
{
	return d_rmul(sqrt(d_abs(z)), d_exp(0.5 * atan2(z.i, z.r)));
}

static inline d_complex_t d_inv(d_complex_t z)
{
	return d_div(d_complex(1.0, 0.0), z);
}

static inline d_complex_t d_neg(d_complex_t z)
{
	return d_complex(-z.r, -z.i);
}

// float

static inline f_complex_t f_complex(double r, double i)
{
	f_complex_t z;

	z.r = (float)r;
	z.i = (float)i;

	return z;
}

static inline f_complex_t f_add(f_complex_t a, f_complex_t b)
{
	return f_complex(a.r + b.r, a.i + b.i);
}

static inline f_complex_t f_sub(f_complex_t a, f_complex_t b)
{
	return f_complex(a.r - b.r, a.i - b.i);
}

static inline f_complex_t f_mul(f_complex_t a, f_complex_t b)
{
	return f_complex((a.r * b.r) - (a.i * b.i), (a.r * b.i) + (a.i * b.r));
}

static inline f_complex_t f_div(f_complex_t a, f_complex_t b)
{
	if ((fabs(b.r) <= 0) && (fabs(b.i) <= 0)) return f_complex(0, 0);
	return f_complex(((a.r * b.r) + (a.i * b.i)) / ((b.r * b.r) + (b.i * b.i)),
	                 ((a.i * b.r) - (a.r * b.i)) / ((b.r * b.r) + (b.i * b.i)));
}

static inline f_complex_t f_rmul(float r, f_complex_t z)
{
	return f_complex(r * z.r, r * z.i);
}

static inline f_complex_t f_add3(f_complex_t a, f_complex_t b, f_complex_t c)
{
	return f_complex(a.r + b.r + c.r, a.i + b.i + c.i);
}

static inline f_complex_t f_add4(f_complex_t a, f_complex_t b, f_complex_t c, f_complex_t d)
{
	return f_complex(a.r + b.r + c.r + d.r, a.i + b.i + c.i + d.i);
}

// cast

static inline d_complex_t d_cast(f_complex_t z)
{
	return d_complex(z.r, z.i);
}

static inline f_complex_t f_cast(d_complex_t z)
{
	return f_complex((float)z.r, (float)z.i);
}

#endif  // _COMPLEX_H_
