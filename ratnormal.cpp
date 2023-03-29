static const char usagetext[] = R"(
Ss{NAME}
    Nm{ratnormal} -- compute rational normal form of a sum of rational
    expressions.

Ss{SYNOPSYS}
    Nm{ratnormal} [Fl{-mndhC}] [Fl{-j} Ar{threads}] Ar{input} [Ar{output}]

Ss{DESCRIPTION}

    Nm{ratnormal} reads an input file, which must contain a sum of
    rational expressions (products of integer powers of polynomials with
    integer coefficients), and brings the sum under a common denominator,
    cancelling the fractions if needed, and preserving the factorization
    of the denominator.
    
    Nm{ratnormal} uses the FLINT library.

Ss{EXAMPLES}

    $ echo '2/3/(x + 2*y^2)^2 - 1/(x + 6*x*y^2)' | ratnormal -
    (
     1/3 *
     1/(6*x*y^2+x) *
     1/(x+2*y^2)^2 *
     (-3*x^2+2*x-12*y^4)
    )

Ss{OPTIONS}
    Fl{-n}         Allow expanding the final numerator to improve performance.
    Fl{-d}         Allow expanding the final denominator.
    Fl{-m}         Take the common monomials out of the brackets at the end.
    Fl{-j} Ar{threads} Use this many worker threads, if possible (default: 1).
    Fl{-h}         Show this help message.
    Fl{-C}         Force colored output even if stderr is not a tty.
    Fl{-V}         Print version information.

Ss{ARGUMENTS}
    Ar{input}      Read the input expression from this file, with "-"
               meaning the standard input.
    Ar{output}     Write the result into this file, with "-" meaning
               the standard output (which is the default).

Ss{AUTHORS}
    Vitaly Magerya <vitaly.magerya@tx97.net>
)";

#include <assert.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <fstream>
#include <ginac/ginac.h>
#include <ginac/parser.h>
#include <iostream>
#include <printf.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static bool COLORS = !!isatty(STDERR_FILENO);

/* Misc
 */

#define SWAP(type, a, b) { type t = a; a = b; b = t; }

static double
timestamp(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

/* Logging
 */

static int log_depth = 0;
static double log_first_timestamp = timestamp();
static double log_last_timestamp = log_first_timestamp;

typedef struct {
    const char *name;
    double time;
} log_func_info;

#define logd(fmt, ...) {\
    double NOW = timestamp(); \
    fprintf(stderr, "[%.4fs +%.4fs]%*s " fmt "\n", NOW - log_first_timestamp, NOW - log_last_timestamp, log_depth*2, "", ## __VA_ARGS__); \
    log_last_timestamp = NOW; \
}

static void
log_func_start(const log_func_info *i)
{
    logd("> %s()", i->name);
    log_depth++;
}

static void
log_func_end(const log_func_info *i)
{
    log_depth--;
    logd("< %s(%.4fs)", i->name, NOW - i->time);
}

#define LOGME \
    __attribute__((cleanup(log_func_end))) log_func_info _log_func_info = {__func__, timestamp()}; \
    log_func_start(&_log_func_info);

/* Rational context
 */

// struct rat_ctx {
fmpz_mpoly_ctx_t ctx;
std::map<std::string, int> variable2index;
std::vector<std::string> variables;
const char **variable_names;
bool factor_numerator = false;
bool factor_denominator = false;
#define nvariables (variables.size())
// }

void
rat_ctx_init()
{
    fmpz_mpoly_ctx_init(ctx, nvariables, ORD_LEX);
}

void
rat_ctx_clear()
{
    fmpz_mpoly_ctx_clear(ctx);
    for (ulong i = 0; i < nvariables; i++) {
        flint_free((void*)variable_names[i]);
        variable_names[i] = NULL;
    }
}

/* Rational
 */

typedef struct {
    fmpq_t numfactor;
    fmpz_mpoly_struct *factors;
    int *powers;
    ulong num;
    ulong alloc;
} rat_struct;

typedef rat_struct rat_t[1];

void
rat_init(rat_t rat)
{
    fmpq_init(rat->numfactor);
    rat->factors = NULL;
    rat->powers = NULL;
    rat->num = 0;
    rat->alloc = 0;
}

void
rat_zero(rat_t rat)
{
    fmpq_zero(rat->numfactor);
    for (ulong i = 0; i < rat->num; i++) {
        fmpz_mpoly_clear(&rat->factors[i], ctx);
    }
    rat->num = 0;
}

bool
rat_is_zero(rat_t rat)
{
    return fmpq_is_zero(rat->numfactor);
}

void
rat_one(rat_t rat)
{
    fmpq_one(rat->numfactor);
    for (ulong i = 0; i < rat->num; i++) {
        fmpz_mpoly_clear(&rat->factors[i], ctx);
    }
    rat->num = 0;
}

void
rat_swap(rat_t a, rat_t b)
{
   rat_struct t = *a;
   *a = *b;
   *b = t;
}

void
rat_clear(rat_t rat)
{
    fmpq_clear(rat->numfactor);
    for (ulong i = 0; i < rat->num; i++) {
        fmpz_mpoly_clear(&rat->factors[i], ctx);
    }
    flint_free(rat->factors);
    flint_free(rat->powers);
    // if paranoid then {
    rat->factors = NULL;
    rat->powers = NULL;
    rat->num = 0;
    rat->alloc = 0;
    // }
}

void
rat_fit_length(rat_t rat, ulong len)
{
    if (rat->alloc < len) {
        len = FLINT_MAX(len, rat->alloc + rat->alloc/2);
        rat->factors = (fmpz_mpoly_struct*)flint_realloc(rat->factors, len*sizeof(rat->factors[0]));
        rat->powers = (int*)flint_realloc(rat->powers, len*sizeof(rat->powers[0]));
        rat->alloc = len;
    }
}

void
rat_set(rat_t rat, const rat_t src)
{
    rat_one(rat);
    fmpq_set(rat->numfactor, src->numfactor);
    rat_fit_length(rat, src->num);
    for (ulong i = 0; i < src->num; i++) {
        fmpz_mpoly_init(&rat->factors[i], ctx);
        fmpz_mpoly_set(&rat->factors[i], &src->factors[i], ctx);
        rat->powers[i] = src->powers[i];
    }
    rat->num = src->num;
}

void
rat_mul_fmpq(rat_t rat, const fmpq_t x, const int power)
{
    assert(!fmpq_is_zero(x));
    fmpq_t xpow;
    fmpq_init(xpow);
    fmpq_pow_si(xpow, x, power);
    fmpq_mul(rat->numfactor, rat->numfactor, xpow);
    fmpq_clear(xpow);
}

void
rat_mul_rat(rat_t rat, const rat_t f)
{
    fmpq_mul(rat->numfactor, rat->numfactor, f->numfactor);
    rat_fit_length(rat, rat->num + f->num);
    for (ulong i = 0; i < f->num; i++) {
        fmpz_mpoly_init(&rat->factors[rat->num + i], ctx);
        fmpz_mpoly_set(&rat->factors[rat->num + i], &f->factors[i], ctx);
        rat->powers[rat->num + i] = f->powers[i];
    }
    rat->num += f->num;
}

void
rat_mul_fmpz(rat_t rat, const fmpz_t x, const int power)
{
    assert(!fmpz_is_zero(x));
    if (power == 0) {
        return;
    } else {
        fmpz_t xpow;
        fmpz_init(xpow);
        if (power > 0) {
            fmpz_pow_ui(xpow, x, power);
            fmpq_mul_fmpz(rat->numfactor, rat->numfactor, xpow);
        } else {
            fmpz_pow_ui(xpow, x, -power);
            fmpq_div_fmpz(rat->numfactor, rat->numfactor, xpow);
        }
        fmpz_clear(xpow);
    }
}

void
rat_mul_fmpz_mpoly(rat_t rat, const fmpz_mpoly_t poly, int power)
{
    if (power == 0) {
        return;
    } else if (fmpz_mpoly_is_fmpz(poly, ctx)) {
        fmpz_t k;
        fmpz_mpoly_get_fmpz(k, poly, ctx);
        assert(!fmpz_is_zero(k));
        rat_mul_fmpz(rat, k, power);
    } else {
        rat_fit_length(rat, rat->num + 1);
        fmpz_mpoly_init(&rat->factors[rat->num], ctx);
        fmpz_mpoly_set(&rat->factors[rat->num], poly, ctx);
        rat->powers[rat->num] = power;
        rat->num++;
    }
}

void
rat_mul_fmpz_mpoly_setx(rat_t rat, fmpz_mpoly_t poly, int power)
{
    assert(!fmpz_mpoly_is_zero(poly, ctx));
    if (power == 0) {
        return;
    } else if (fmpz_mpoly_is_fmpz(poly, ctx)) {
        fmpz_t k;
        fmpz_init(k);
        fmpz_mpoly_get_fmpz(k, poly, ctx);
        rat_mul_fmpz(rat, k, power);
        fmpz_mpoly_zero(poly, ctx);
    } else {
        rat_fit_length(rat, rat->num + 1);
        fmpz_mpoly_init(&rat->factors[rat->num], ctx);
        fmpz_mpoly_swap(&rat->factors[rat->num], poly, ctx);
        rat->powers[rat->num] = power;
        rat->num++;
    }
}

void
rat_fprint(FILE *f, const rat_t rat)
{
    if ((rat->num == 0) || !fmpq_is_one(rat->numfactor)) {
        fprintf(f, " ");
        fmpq_fprint(f, rat->numfactor);
        if (0 != rat->num) {
            fprintf(f, " *\n");
        }
    }
    for (ulong i = 0; i < rat->num; i++) {
        int power = rat->powers[i];
        if (power >= 0) {
            fprintf(f, " (");
        } else {
            fprintf(f, " 1/(");
            power = -power;
        }
        fmpz_mpoly_fprint_pretty(f, &rat->factors[i], variable_names, ctx);
        if (power != 1) {
            fprintf(f, ")^%d", power);
        } else {
            fprintf(f, ")");
        }
        if (i != rat->num-1) {
            fprintf(f, " *\n");
        }
    }
}

void
rat_fprint_short(FILE *f, const rat_t rat)
{
    fprintf(f, "R{");
    for (ulong i = 0; i < rat->num; i++) {
        slong nterms = fmpz_mpoly_length(&rat->factors[i], ctx);
        int power = rat->powers[i];
        if (i == 0) {
            if (power < 0) {
                fprintf(f, "1/");
                power = -power;
            }
        } else {
            if (power < 0) {
                fprintf(f, "/");
                power = -power;
            } else {
                fprintf(f, " ");
            }
        }
        if (power == 1) {
            fprintf(f, "%ldt", nterms);
        } else {
            fprintf(f, "%ldt^%d", nterms, power);
        }
    }
    fprintf(f, "}");
}

void
rat_swapoff_fmpz_factor(rat_t rat, ulong i)
{
    assert(i < rat->num);
    assert(fmpz_mpoly_is_fmpz(&rat->factors[i], ctx));
    rat_mul_fmpz_mpoly_setx(rat, &rat->factors[i], rat->powers[i]);
    if (i < rat->num-1) {
        fmpz_mpoly_swap(&rat->factors[i], &rat->factors[rat->num-1], ctx);
        SWAP(int, rat->powers[i], rat->powers[rat->num-1]);
    }
    fmpz_mpoly_clear(&rat->factors[rat->num-1], ctx);
    rat->num--;
}

slong
rat_max_length(const rat_t rat)
{
    slong maxlen = 0;
    for (ulong i = 0; i < rat->num; i++) {
        int pow = rat->powers[i];
        if (pow < 0) pow = -pow;
        maxlen = FLINT_MAX(maxlen, pow*fmpz_mpoly_length(&rat->factors[i], ctx));
    }
    return maxlen;
}

void
rat_sort(rat_t rat)
{
    for (ulong i = 0; i < rat->num; i++) {
        slong nterms1 = fmpz_mpoly_length(&rat->factors[i], ctx);
        for (ulong j = i+1; j < rat->num; j++) {
            slong nterms2 = fmpz_mpoly_length(&rat->factors[j], ctx);
            if (nterms2 > nterms1) {
                fmpz_mpoly_swap(&rat->factors[i], &rat->factors[j], ctx);
                SWAP(int, rat->powers[i], rat->powers[j]);
                nterms1 = nterms2;
            }
        }
    }
}

void
rat_reverse(rat_t rat)
{
    for (ulong i = 0; i < rat->num/2; i++) {
        ulong j = rat->num - 1 - i;
        assert(i < j);
        fmpz_mpoly_swap(&rat->factors[i], &rat->factors[j], ctx);
        SWAP(int, rat->powers[i], rat->powers[j]);
    }
}

void
rat_cofactorize(rat_t rat)
{
    LOGME;
    logd("Cofactorizing %r", rat);
    fmpz_mpoly_t gcd;
    fmpz_mpoly_init(gcd, ctx);
    for (ulong i = 0; i < rat->num; i++) {
        #define A (&rat->factors[i])
        #define Apower rat->powers[i]
        //if (Apower >= 0) continue;
        for (ulong j = i + 1; j < rat->num; j++) {
            #define B (&rat->factors[j])
            #define Bpower rat->powers[j]
            //if (Bpower >= 0) continue;
            fmpz_mpoly_gcd(gcd, A, B, ctx);
            if (fmpz_mpoly_is_one(gcd, ctx)) {
                continue;
            }
            logd("Split factor %p^%d from %p^%d %p^%d", gcd, Apower + Bpower, A, Apower, B, Bpower);
            if (fmpz_mpoly_divides(A, A, gcd, ctx) != 1) { assert(0); }
            if (fmpz_mpoly_divides(B, B, gcd, ctx) != 1) { assert(0); }
            rat_mul_fmpz_mpoly_setx(rat, gcd, Apower + Bpower);
            if (fmpz_mpoly_is_fmpz(B, ctx)) {
                rat_swapoff_fmpz_factor(rat, j);
                j--;
            }
            if (fmpz_mpoly_is_fmpz(A, ctx)) {
                rat_swapoff_fmpz_factor(rat, i);
                i--;
                break;
            }
            #undef B
            #undef Bpower
        }
        #undef A
        #undef Apower
    }
    fmpz_mpoly_clear(gcd, ctx);
}

void
poly_reduce_mul(fmpz_mpoly_struct *factors, ulong n)
{
    if (n == 0) { return; }
    for (; n > 1; n--) {
        // Find two shortest polynomials, multiply them.
        slong length1 = LONG_MAX, idx1 = -1;
        slong length2 = LONG_MAX, idx2 = -1;
        for (ulong i = 0; i < n; i++) {
            slong length = fmpz_mpoly_length(&factors[i], ctx);
            if (length < length1) {
                length2 = length1; idx2 = idx1;
                length1 = length; idx1 = i;
            } else if (length < length2) {
                length2 = length; idx2 = i;
            }
        }
        assert((idx1 >= 0) && (idx2 >= 0) && (idx1 != idx2));
        logd("Multiplying %p*%p", &factors[idx1], &factors[idx2]);
        fmpz_mpoly_mul(&factors[idx1], &factors[idx1], &factors[idx2], ctx);
        fmpz_mpoly_clear(&factors[idx2], ctx);
        factors[idx2] = factors[n-1];
    }
}

void
rat_expand_numerator(fmpz_mpoly_t poly, const rat_t rat)
{
    LOGME;
    logd("Expanding the numerator of %r", rat);
    std::vector<fmpz_mpoly_struct> factors;
    for (ulong i = 0; i < rat->num; i++) {
        if (rat->powers[i] <= 0) continue;
        fmpz_mpoly_t tmp;
        fmpz_mpoly_init(tmp, ctx);
        if (rat->powers[i] > 1) {
            logd("Expanding %p^%d", &rat->factors[i], rat->powers[i]);
            fmpz_mpoly_pow_ui(tmp, &rat->factors[i], rat->powers[i], ctx);
        } else {
            fmpz_mpoly_set(tmp, &rat->factors[i], ctx);
        }
        factors.push_back(*tmp);
    }
    slong n = factors.size();
    if (n == 0) {
        fmpz_mpoly_one(poly, ctx);
    } else {
        poly_reduce_mul(&factors[0], factors.size());
        fmpz_mpoly_clear(poly, ctx);
        fmpz_mpoly_set(poly, &factors[0], ctx);
    }
    if (!fmpz_is_one(fmpq_numref(rat->numfactor))) {
        fmpz_mpoly_scalar_mul_fmpz(poly, poly, fmpq_numref(rat->numfactor), ctx);
    }
}

void
rat_expand_denominator(fmpz_mpoly_t poly, const rat_t rat)
{
    LOGME;
    logd("Expanding the denominator of %r", rat);
    std::vector<fmpz_mpoly_struct> factors;
    for (ulong i = 0; i < rat->num; i++) {
        if (rat->powers[i] >= 0) continue;
        fmpz_mpoly_t tmp;
        fmpz_mpoly_init(tmp, ctx);
        if (rat->powers[i] < -1) {
            logd("Expanding %p^%d", &rat->factors[i], -rat->powers[i]);
            fmpz_mpoly_pow_ui(tmp, &rat->factors[i], -rat->powers[i], ctx);
        } else {
            fmpz_mpoly_set(tmp, &rat->factors[i], ctx);
        }
        factors.push_back(*tmp);
    }
    slong n = factors.size();
    if (n == 0) {
        fmpz_mpoly_one(poly, ctx);
    } else {
        poly_reduce_mul(&factors[0], factors.size());
        fmpz_mpoly_clear(poly, ctx);
        fmpz_mpoly_set(poly, &factors[0], ctx);
    }
    if (!fmpz_is_one(fmpq_denref(rat->numfactor))) {
        fmpz_mpoly_scalar_mul_fmpz(poly, poly, fmpq_denref(rat->numfactor), ctx);
    }
}

void
rat_add_setx(rat_t res, rat_t rat1, rat_t rat2)
{
    LOGME;
    logd("Adding %r and %r", rat1, rat2);
    if (rat_is_zero(rat1)) {
        rat_clear(res);
        *res = *rat2;
        return;
    }
    if (rat_is_zero(rat2)) {
        rat_clear(res);
        *res = *rat1;
        return;
    }
    rat_one(res);
    // First, take out common factors.
    fmpz_mpoly_t gcd, aprime, bprime;
    fmpz_mpoly_init(gcd, ctx);
    fmpz_mpoly_init(aprime, ctx);
    fmpz_mpoly_init(bprime, ctx);
    #define A (&rat1->factors[i])
    #define Apower rat1->powers[i]
    #define B (&rat2->factors[j])
    #define Bpower rat2->powers[j]
    for (ulong i = 0; i < rat1->num; i++) {
        if (Apower >= 0 && factor_numerator) continue;
        for (ulong j = 0; j < rat2->num; j++) {
            if (Bpower >= 0 && factor_numerator) continue;
            fmpz_mpoly_gcd_cofactors(gcd, aprime, bprime, A, B, ctx);
            if (fmpz_mpoly_is_one(gcd, ctx)) {
                continue;
            }
            int power = FLINT_MAX(Apower, Bpower);
            //logd("Common factor: %p^%d", gcd, power);
            fmpz_mpoly_swap(A, aprime, ctx);
            fmpz_mpoly_swap(B, bprime, ctx);
            rat_mul_fmpz_mpoly(rat1, gcd, Apower - power);
            rat_mul_fmpz_mpoly(rat2, gcd, Bpower - power);
            rat_mul_fmpz_mpoly_setx(res, gcd, power);
            if (fmpz_mpoly_is_fmpz(B, ctx)) {
                rat_swapoff_fmpz_factor(rat2, j);
                j--;
            }
            if (fmpz_mpoly_is_fmpz(A, ctx)) {
                rat_swapoff_fmpz_factor(rat1, i);
                i--;
                break;
            }
        }
    }
    fmpz_mpoly_clear(gcd, ctx);
    fmpz_mpoly_clear(aprime, ctx);
    fmpz_mpoly_clear(bprime, ctx);
    // Don't forget about the numeric factors.
    fmpz_t g;
    fmpz_init(g);
    fmpz_gcd(g, fmpq_denref(rat1->numfactor), fmpq_denref(rat2->numfactor));
    assert(!fmpz_is_zero(g));
    fmpq_mul_fmpz(rat1->numfactor, rat1->numfactor, g);
    fmpq_mul_fmpz(rat2->numfactor, rat2->numfactor, g);
    rat_mul_fmpz(res, g, -1);
    fmpz_clear(g);
    logd("Common factor: %r", res);
    logd("Remaining: %r and %r", rat1, rat2);
    // Then, expand the combined numerator.
    fmpz_mpoly_t anum, aden, bnum, bden;
    fmpz_mpoly_init(anum, ctx);
    fmpz_mpoly_init(bnum, ctx);
    fmpz_mpoly_init(aden, ctx);
    fmpz_mpoly_init(bden, ctx);
    rat_expand_numerator(anum, rat1);
    assert(!fmpz_mpoly_is_zero(anum, ctx));
    rat_expand_denominator(aden, rat1);
    assert(!fmpz_mpoly_is_zero(aden, ctx));
    rat_expand_numerator(bnum, rat2);
    assert(!fmpz_mpoly_is_zero(bnum, ctx));
    rat_expand_denominator(bden, rat2);
    assert(!fmpz_mpoly_is_zero(bden, ctx));
    logd("Computing %p*%p + %p*%p", anum, bden, bnum, aden);
    fmpz_mpoly_mul(anum, anum, bden, ctx);
    fmpz_mpoly_mul(bnum, bnum, aden, ctx);
    fmpz_mpoly_add(anum, anum, bnum, ctx);
    if (fmpz_mpoly_is_zero(anum, ctx)) {
        rat_zero(res);
        fmpz_mpoly_clear(anum, ctx);
        fmpz_mpoly_clear(bnum, ctx);
        rat_clear(rat1);
        rat_clear(rat2);
        fmpz_mpoly_clear(aden, ctx);
        fmpz_mpoly_clear(bden, ctx);
        return;
    }
    rat_mul_fmpz_mpoly_setx(res, anum, 1);
    fmpz_mpoly_clear(anum, ctx);
    fmpz_mpoly_clear(bnum, ctx);
    // Finally, multiply by the combined denominator.
    if (!factor_denominator) {
        for (ulong i = 0; i < rat1->num; i++) {
            if (Apower >= 0) continue;
            rat_mul_fmpz_mpoly_setx(res, A, Apower);
        }
        for (ulong j = 0; j < rat2->num; j++) {
            if (Bpower >= 0) continue;
            rat_mul_fmpz_mpoly_setx(res, B, Bpower);
        }
        rat_mul_fmpz(res, fmpq_denref(rat1->numfactor), -1);
        rat_mul_fmpz(res, fmpq_denref(rat2->numfactor), -1);
    } else {
        rat_mul_fmpz_mpoly_setx(res, aden, -1);
        rat_mul_fmpz_mpoly_setx(res, bden, -1);
    }
    fmpz_mpoly_clear(aden, ctx);
    fmpz_mpoly_clear(bden, ctx);
    #undef A
    #undef Apower
    #undef B
    #undef Bpower
    rat_sort(res);
    rat_cofactorize(res);
    rat_sort(res);
    logd("Result is %r", res);
}

void
rat_sqrt(rat_t rat)
{
    LOGME;
    logd("Sqrt of %r", rat);
    rat_cofactorize(rat);
    bool minus = false;
    if (fmpq_sgn(rat->numfactor) < 0) {
        minus = !minus;
        fmpq_neg(rat->numfactor, rat->numfactor);
    }
    fmpz_sqrt(fmpq_numref(rat->numfactor), fmpq_numref(rat->numfactor));
    fmpz_sqrt(fmpq_denref(rat->numfactor), fmpq_denref(rat->numfactor));
    fmpz_mpoly_t root;
    fmpz_mpoly_init(root, ctx);
    for (ulong i = 0; i < rat->num; i++) {
        int n = rat->powers[i];
        if ((n & 1) == 0) {
            rat->powers[i] = n / 2;
        } else {
            logd("Sqrt of %p", &rat->factors[i]);
            if (!fmpz_mpoly_sqrt(root, &rat->factors[i], ctx)) {
                minus = !minus;
                fmpz_mpoly_neg(&rat->factors[i], &rat->factors[i], ctx);
                if (!fmpz_mpoly_sqrt(root, &rat->factors[i], ctx)) {
                    fprintf(stderr, "ratnormal: sqrt() of a non-square\n");
                    exit(1);
                }
            }
            fmpz_mpoly_swap(root, &rat->factors[i], ctx);
        }
    }
    fmpz_mpoly_clear(root, ctx);
    if (minus) {
        fprintf(stderr, "ratnormal: sqrt() of a negative number\n");
        exit(1);
    }
}

void
rat_take_out_monomials(rat_t rat)
{
    fmpz_mpoly_t m;
    fmpz_mpoly_init(m, ctx);
    for (ulong i = 0; i < rat->num; i++) {
        if (fmpz_mpoly_length(&rat->factors[i], ctx) < 2) continue;
        fmpz_mpoly_term_content(m, &rat->factors[i], ctx);
        if (fmpz_mpoly_is_one(m, ctx)) continue;
        if (fmpz_mpoly_divides(&rat->factors[i], &rat->factors[i], m, ctx) != 1) { assert(0); }
        if (rat->powers[i] > 0) {
            if (fmpz_mpoly_pow_ui(m, m, rat->powers[i], ctx) != 1) { assert(0); };
            rat_mul_fmpz_mpoly_setx(rat, m, 1);
        } else {
            if (fmpz_mpoly_pow_ui(m, m, -rat->powers[i], ctx) != 1) { assert(0); };
            rat_mul_fmpz_mpoly_setx(rat, m, -1);
        }
    }
    fmpz_mpoly_clear(m, ctx);
}

void
rat_reduce_add(rat_struct *rats, ulong n)
{
    LOGME;
    for (; n > 1; n--) {
        // Find two shortest rationals, add them.
        slong length1 = LONG_MAX, idx1 = -1;
        slong length2 = LONG_MAX, idx2 = -1;
        for (ulong i = 0; i < n; i++) {
            slong length = rat_max_length(&rats[i]);
            if (length < length1) {
                length2 = length1; idx2 = idx1;
                length1 = length; idx1 = i;
                continue;
            }
            if (length < length2) {
                length2 = length; idx2 = i;
                continue;
            }
        }
        assert((idx1 >= 0) && (idx2 >= 0) && (idx1 != idx2));
        logd("Adding %r and %r (todo: %lu)", &rats[idx1], &rats[idx2], n-2);
        rat_t rat;
        rat_init(rat);
        rat_add_setx(rat, &rats[idx1], &rats[idx2]);
        rat_swap(&rats[idx1], rat);
        rat_swap(&rats[idx2], &rats[n-1]);
    }
}

/* GiNaC conversion
 */

template <typename F> void
factor_iter(const GiNaC::ex &e, F yield)
{
    if (GiNaC::is_a<GiNaC::mul>(e)) {
        for (const auto &f : e) {
            if (GiNaC::is_a<GiNaC::power>(f)) {
                yield(f.op(0), GiNaC::ex_to<GiNaC::numeric>(f.op(1)).to_int());
            } else {
                yield(f, 1);
            }
        }
    } else {
        if (GiNaC::is_a<GiNaC::power>(e)) {
            yield(e.op(0), GiNaC::ex_to<GiNaC::numeric>(e.op(1)).to_int());
        } else {
            yield(e, 1);
        }
    }
}

template <typename F> void
term_iter(const GiNaC::ex &e, F yield)
{
    if (GiNaC::is_a<GiNaC::add>(e)) {
        for (const auto &t : e) {
            yield(t);
        }
    } else {
        yield(e);
    }
}

void
fmpz_of_ginac(fmpz_t x, const GiNaC::numeric &num)
{
    assert(num.is_integer());
    assert((LONG_MIN <= num) && (num <= LONG_MAX));
    fmpz_set_si(x, num.to_long());
}

void
fmpq_of_ginac(fmpq_t x, const GiNaC::numeric &n)
{
    assert(n.is_rational());
    const GiNaC::numeric num = n.numer();
    const GiNaC::numeric den = n.denom();
    assert((LONG_MIN <= num) && (num <= LONG_MAX));
    assert((LONG_MIN <= den) && (den <= LONG_MAX));
    fmpq_set_si(x, num.to_long(), den.to_long());
}

void rat_of_ginac(rat_t rat, const GiNaC::ex &expr);

void
oldrat_of_ginac(rat_t rat, const GiNaC::ex &expr)
{
    LOGME;
    rat_one(rat);
    std::vector<ulong> exp(nvariables);
    factor_iter(expr, [&](const GiNaC::ex &polyfactor, int pfpower) {
        if (GiNaC::is_a<GiNaC::numeric>(polyfactor)) {
            GiNaC::numeric npf = GiNaC::ex_to<GiNaC::numeric>(polyfactor);
            assert(npf.is_rational());
            fmpq_t n;
            fmpq_init(n);
            fmpq_of_ginac(n, npf);
            rat_mul_fmpq(rat, n, pfpower);
            fmpq_clear(n);
        } else {
            fmpz_mpoly_t poly;
            fmpz_mpoly_init(poly, ctx);
            term_iter(polyfactor.expand(), [&](const GiNaC::ex &term) {
                fmpz_t coef;
                fmpz_init_set_ui(coef, 1);
                for (ulong i = 0; i < nvariables; i++) {
                    exp[i] = 0;
                }
                factor_iter(term, [&](const GiNaC::ex &f, int tfpower) {
                    assert(tfpower >= 0);
                    if (GiNaC::is_a<GiNaC::numeric>(f)) {
                        GiNaC::numeric ntf = GiNaC::ex_to<GiNaC::numeric>(f);
                        assert(ntf.is_integer());
                        fmpz_t npow;
                        fmpz_init(npow);
                        fmpz_of_ginac(npow, ntf);
                        fmpz_pow_ui(npow, npow, tfpower);
                        fmpz_mul(coef, coef, npow);
                        fmpz_clear(npow);
                    } else if (GiNaC::is_a<GiNaC::symbol>(f)) {
                        GiNaC::symbol sym = GiNaC::ex_to<GiNaC::symbol>(f);
                        int varidx = variable2index[sym.get_name()];
                        exp[varidx] += tfpower;
                    } else {
                        assert(!"A term has a factor which is neither a number nor a symbol");
                    }
                });
                fmpz_mpoly_push_term_fmpz_ui(poly, coef, &exp[0], ctx);
                fmpz_clear(coef);
            });
            fmpz_mpoly_sort_terms(poly, ctx);
            fmpz_mpoly_combine_like_terms(poly, ctx);
            rat_mul_fmpz_mpoly_setx(rat, poly, pfpower);
            fmpz_mpoly_clear(poly, ctx);
        }
    });
}

/* Main
 */

int
print_ptr_arginfo(const struct printf_info *info, size_t n, int *argtypes)
{
    if (n > 0) argtypes[0] = PA_POINTER;
    return 1;
}

int
print_fmpz(FILE *f, const struct printf_info *info, const void *const *args)
{
    fmpz_fprint(f, **((const fmpz_t **)(args[0])));
    return 1;
}

int
print_poly(FILE *f, const struct printf_info *info, const void *const *args)
{
    fmpz_mpoly_fprint_pretty(f, **((const fmpz_mpoly_t **)(args[0])), variable_names, ctx);
    return 1;
}

int
print_poly_short(FILE *f, const struct printf_info *info, const void *const *args)
{
    slong nterms = fmpz_mpoly_length(**((const fmpz_mpoly_t **)(args[0])), ctx);
    fprintf(f, "P{%ldt}", nterms);
    return 1;
}

int
print_rat(FILE *f, const struct printf_info *info, const void *const *args)
{
    rat_fprint(f, **((const rat_t **)(args[0])));
    return 1;
}

int
print_rat_short(FILE *f, const struct printf_info *info, const void *const *args)
{
    rat_fprint_short(f, **((const rat_t **)(args[0])));
    return 1;
}

// GiNaC will crash if these functions are not aligned...
static __attribute__((aligned(8))) GiNaC::ex
sqrt_reader(const GiNaC::exvector& ev)
{ return GiNaC::sqrt(ev[0]); }

static __attribute__((aligned(8))) GiNaC::ex
power_reader(const GiNaC::exvector& ev)
{ return GiNaC::power(ev[0], ev[1]); }

GiNaC::ex
load_input(const char *filename)
{
    LOGME;
    GiNaC::prototype_table proto = GiNaC::prototype_table();
    proto[std::make_pair("Sqrt", 1)] = sqrt_reader;
    proto[std::make_pair("sqrt", 1)] = sqrt_reader;
    proto[std::make_pair("Power", 2)] = power_reader;
    proto[std::make_pair("pow", 2)] = power_reader;
    GiNaC::parser reader(GiNaC::symtab(), false, proto);
    GiNaC::ex expr;
    if (strcmp(filename, "-") != 0) {
        std::ifstream ifs(filename);
        expr = reader(ifs);
    } else {
        expr = reader(std::cin);
    }
    // List variable names, init context.
    auto &name2sym = reader.get_syms();
    variable_names = (const char**)flint_calloc(name2sym.size(), sizeof(const char*));
    for (auto &&kv : name2sym) {
        int index = nvariables;
        variables.push_back(kv.first);
        variable_names[index] = (const char*)flint_calloc(variables[index].length() + 1, 1);
        strcpy((char*)variable_names[index], variables[index].c_str());
        variable2index[kv.first] = index;
        logd("Variable %d is %s", index+1, variable_names[index]);
    }
    rat_ctx_init();
    return expr;
}

void
save_output(const char *filename, rat_t rat)
{
    LOGME;
    logd("Saving %r to %s", rat, filename);
    FILE *f;
    if (strcmp(filename, "-") != 0) {
        f = fopen(filename, "wb");
        if (f == NULL) {
            fprintf(stderr, "ratnormal: can't open file %s\n", filename);
            exit(1);
        }
    } else {
        f = stdout;
    }
    fprintf(f, "(\n");
    rat_fprint(f, rat);
    fprintf(f, "\n)\n");
    if (strcmp(filename, "-") != 0) {
        fclose(f);
    } else {
        fflush(f);
    }
}

void
usage(FILE *f)
{
    const char *p = strchr(usagetext, '\n') + 1;
    for (;;) {
        const char *l = strchr(p + 2, '{');
        if (l == NULL) break;
        const char *r = strchr(l, '}');
        if (r == NULL) break;
        const char *a = "", *b = "\033[0m";
        if (l[-2] == 'S' && l[-1] == 's') { a = "\033[1m"; }
        if (l[-2] == 'N' && l[-1] == 'm') { a = "\033[1;35m"; }
        if (l[-2] == 'F' && l[-1] == 'l') { a = "\033[33m"; }
        if (l[-2] == 'C' && l[-1] == 'm') { a = "\033[1m"; }
        if (l[-2] == 'A' && l[-1] == 'r') { a = "\033[32m"; }
        if (l[-2] == 'D' && l[-1] == 'l') { a = "\033[36m"; }
        if (l[-2] == 'Q' && l[-1] == 'l') { a = "\033[36m"; }
        fwrite(p, l - p - 2, 1, f);
        if (COLORS) fputs(a, f);
        fwrite(l + 1, r - l - 1, 1, f);
        if (COLORS) fputs(b, f);
        p = r + 1;
    }
    fputs(p, f);
}

void
rat_pow(rat_t rat, int pow)
{
    fmpq_pow_si(rat->numfactor, rat->numfactor, pow);
    for (ulong i = 0; i < rat->num; i++) {
        rat->powers[i] *= pow;
    }
}

void
rat_of_ginac(rat_t rat, const GiNaC::ex &expr)
{
    //LOGME;
    if (GiNaC::is_a<GiNaC::add>(expr)) {
        size_t n = expr.nops();
        assert(n > 0);
        std::vector<rat_struct> rats;
        rats.resize(n);
        rats[0] = *rat;
        for (size_t i = 0; i < n; i++) {
            if (i) rat_init(&rats[i]);
            rat_of_ginac(&rats[i], expr.op(i));
        }
        if (n > 1) rat_reduce_add(&rats[0], n);
        *rat = rats[0];
        return;
    }
    if (GiNaC::is_a<GiNaC::mul>(expr)) {
        size_t n = expr.nops();
        assert(n > 0);
        rat_t f;
        rat_init(f);
        rat_one(rat);
        for (size_t i = 0; i < n; i++) {
            rat_of_ginac(f, expr.op(i));
            rat_mul_rat(rat, f);
        }
        return;
    }
    if (GiNaC::is_a<GiNaC::numeric>(expr)) {
        GiNaC::numeric num = GiNaC::ex_to<GiNaC::numeric>(expr);
        rat_one(rat);
        fmpq_of_ginac(rat->numfactor, num);
        return;
    }
    if (GiNaC::is_a<GiNaC::symbol>(expr)) {
        fmpz_t coef;
        fmpz_init_set_ui(coef, 1);
        std::vector<ulong> exp(nvariables, 0);
        GiNaC::symbol sym = GiNaC::ex_to<GiNaC::symbol>(expr);
        int varidx = variable2index[sym.get_name()];
        exp[varidx] = 1;
        fmpz_mpoly_t poly;
        fmpz_mpoly_init(poly, ctx);
        fmpz_mpoly_push_term_fmpz_ui(poly, coef, &exp[0], ctx);
        fmpz_clear(coef);
        rat_one(rat);
        rat_mul_fmpz_mpoly_setx(rat, poly, 1);
        fmpz_mpoly_clear(poly, ctx);
        return;
    }
    if (GiNaC::is_a<GiNaC::power>(expr)) {
        GiNaC::numeric pow = GiNaC::ex_to<GiNaC::numeric>(expr.op(1));
        if (pow.is_integer()) {
            rat_of_ginac(rat, expr.op(0));
            rat_pow(rat, pow.to_int());
        } else {
            GiNaC::numeric num = pow.numer();
            GiNaC::numeric den = pow.denom();
            if (den == 2) {
                rat_of_ginac(rat, expr.op(0));
                rat_pow(rat, num.to_int());
                rat_sqrt(rat);
            } else {
                fprintf(stderr, "ratnormal: unsupported fractional exponent\n");
                exit(1);
            }
        }
        return;
    }
    fprintf(stderr, "ratnormal: unsupported expression type\n");
    std::cerr << expr;
    exit(1);
}

int
main(int argc, char *argv[])
{
    register_printf_function('F', print_fmpz, print_ptr_arginfo);
    register_printf_function('P', print_poly, print_ptr_arginfo);
    register_printf_function('p', print_poly_short, print_ptr_arginfo);
    register_printf_function('R', print_rat, print_ptr_arginfo);
    register_printf_function('r', print_rat_short, print_ptr_arginfo);
    int nthreads = 1;
    bool take_out_monomials = false;
    const char *inputfile = "-";
    const char *outputfile = "-";
    for (int opt; (opt = getopt(argc, (char*const*)argv, "j:hCVndmS")) != -1;) {
        switch (opt) {
        case 'j': nthreads = atoi(optarg); break;
        case 'h': usage(stdout); return 0;
        case 'V': printf("%s", VERSION); return 0;
        case 'C': COLORS = true; break;
        case 'n': factor_numerator = true; break;
        case 'd': factor_denominator = true; break;
        case 'm': take_out_monomials = true; break;
        default: return 1;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc == 1) {
        inputfile = argv[0];
    } else if (argc == 2) {
        inputfile = argv[0];
        outputfile = argv[1];
    } else {
        fprintf(stderr, "ratnormal: bad invocation (use -h to see usage)\n");
        return 1;
    }
    flint_set_num_threads(nthreads);
    GiNaC::ex expr = load_input(inputfile);
    rat_t rat;
    rat_init(rat);
    rat_of_ginac(rat, expr);
    rat_sort(rat);
    rat_reverse(rat);
    if (take_out_monomials) {
        rat_take_out_monomials(rat);
    }
    save_output(outputfile, rat);
    rat_clear(rat);
    rat_ctx_clear();

    return 0;
}
