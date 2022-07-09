#include "s21_math.h"

long int s21_abs(int x) {
    if (x < 0)
        x *= -1;
    return x;
}
long double s21_fabs(double x) {
    if (x < 0)
        x *= -1;
    return x;
}
// PI > x > PI
double reader(double x) {
    while (x > s21_PI || x < -s21_PI) {
    x += x > s21_PI ? -2 * s21_PI : 2 * s21_PI;
    }
    return x;
}
// sin
long double s21_sin(double x) {
x = reader(x);
    long double rez = x, y = x;
    long double i = 1.;
    while (s21_fabs(rez) > s21_EPS) {
        rez = -1 * rez * x * x / (2 * i * (2 * i + 1));
        i += 1.;
        y += rez;
    }
    return y;
}
// cos
long double s21_cos(double x) {
    x = reader(x);
    return s21_sin((s21_PI / 2.0) - x);
}
// tan
long double s21_tan(double x) {
    x = reader(x);
    return s21_sin(x) / s21_cos(x);
}
// arcsin
long double s21_asin(double x) {
    long double rez = x, y = x;
    long double i = 1;
    while (s21_fabs(rez) > s21_EPS) {
        if (x < -1 || x > 1) {
            y = s21_NAN;
            break;
        }
        if (x == 1 || x == -1) {
            y = s21_PI / 2 * x;
            break;
        }
        rez *=
            x * x * (2 * i - 1) * (2 * i - 1) / ((2 * i + 1) * 2 * i);
        i += 1;
        y += rez;
    }
    return y;
}
// arccos
long double s21_acos(double x) {
    if (x <= 1 && x >= -1) {
        x = s21_PI / 2. - s21_asin(x);
    } else {
        x = s21_NAN;
    }
    return x;
}
// arctan
long double s21_atan(double x) {
    long double y = 0;
    y = (x < 1 && x > -1) ? checker(x) : y;
    y = x == 1 ? s21_PI / 4 : x == -1 ? s21_NPI / 4
                          : x == 0    ? 0
                                      : y;
    y = x > 1 ? s21_PI / 2 - checker(1 / x) : x < -1 ? s21_NPI / 2 - checker(1 / x)
                                                     : y;
    return y;
}
// 1 > arctan > -1
long double checker(double x) {
    long double rez = x, y = x, i = 1;
    while (s21_fabs(rez) > s21_EPS) {
        rez = -1 * rez * x * x * (2 * i - 1) / (2 * i + 1);
        i += 1;
        y += rez;
    }
    return y;
}
// nearest int larger or equal
long double s21_ceil(double x) {
    int count;
    double y = x < 0 ? -x : x;
    for (count = 0; count < y; count++) {
    }
    return x < 0 ? -count + 1 : count;
}
// nearest int less or equal
long double s21_floor(double x) {
    int count;
    double y = x < 0 ? -x : x;
    for (count = 0; count < y; count++) {
    }
    return x < 0 ? -count : x == 0 ? 0
                                   : count - 1;
}
// e raised to xth power
long double s21_exp(double x) {
    long double rez = 1, y = 1;
    long double i = 1;
    int flag = 0;
    if (x < 0) {
        x *= -1;
        flag = 1;
    }
    while (s21_fabs(rez) > s21_EPS) {
    rez *= x / i;
        i++;
        y += rez;
        if (y > s21_MAX_double) {
            y = s21_INF;
            break;
        }
    }
    y = flag == 1 ? y > s21_MAX_double ? 0 : 1. / y : y;
    return y = y > s21_MAX_double ? s21_INF : y;
}
// remainder x / y
long double s21_fmod(double x, double y) {
    long double rez = s21_fabs(x);
    y = s21_fabs(y);
    int i = 1;
    while (rez > y) {
        rez = s21_fabs(x) - y * i;
        i++;
    }
    return x < 0 ? -rez : rez;
}

// logarithm
long double s21_log(double x) {
    long double y;
    return y = (x > 0 && x < 2) ? log_0_2(x) : log_other(x);
}

// logarithm 02
long double log_0_2(double x) {
    x--;
    long double rez = x, y = x;
    long double cnt = 2;
    while (s21_fabs(rez) > s21_EPS) {
        rez *= -x * (cnt - 1) / cnt;
        cnt += 1;
        y += rez;
    }
    return y;
}

// logarithm other
long double log_other(double x) {
    struct special _special;
    translate(x, &_special);
    x = _special.mantissa * _special.sign / 10;
    long double rez;
    rez = x < 0    ? s21_N_NAN
          : x == 0 ? s21_NINF
                   : log_0_2(x) + (_special.e + 1) * s21_ln10;
    return rez;
}

// translate
void translate(double x, struct special *_special) {
    long double i = 1;
    int es = 0;
    _special->sign = x < 0 ? -1 : 1;
    x *= _special->sign;
    if (x >= 10) {
        while (x >= 10) {
            x /= 10.;
            i *= 10;
            es++;
        }
    } else if (x < 1 && x > 0) {
        while (x < 1) {
            x *= 10;
            i /= 10;
            es--;
        }
    }
    _special->mantissa = x;
    _special->pow = i;
    _special->e = es;
}

// base raised to the power of exp
long double s21_pow(double base, double exp) {
    long double rez;
    rez = s21_exp(exp * s21_log(base));
    return rez;
}

// square root of x
long double s21_sqrt(double x) {
    long double rez = 4, y = 0;
    while (s21_fabs(rez - y) > s21_EPS) {
        if (x < 0) {
            rez = s21_N_NAN;
            break;
        }
        y = rez;
        rez = (y + x / y) / 2;
    }
    return rez;
}
